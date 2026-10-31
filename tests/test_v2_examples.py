"""V2 teaching files: compare independently generated CSV to production TMM."""
from pathlib import Path
import numpy as np
from multilayer_forward import coherent_stack

ROOT=Path(__file__).resolve().parents[1]/"examples"/"v2"

def test_v2_synthetic_examples_match_forward_model():
    for filename,d,absorbing in [
        ("transparent_film_120nm.csv",120.,False),
        ("transparent_film_240nm.csv",240.,False),
        ("absorbing_film_180nm.csv",180.,True),
    ]:
        a=np.genfromtxt(ROOT/filename,delimiter=",",names=True)
        w=a["wavelength_nm"]
        assert len(w)==81
        assert np.all(np.diff(w)>0)
        if absorbing:
            n=1.73+0.12*(550/w-1)
            k=0.07*np.exp(-0.5*((w-610)/60)**2)
            index=n+1j*k
        else:
            index=1.65
        R,T=coherent_stack(w,[(index,d)],substrate_n=1.46)
        assert np.max(np.abs(R-a["R"]))<1e-7,(filename,"R")
        assert np.max(np.abs(T-a["T"]))<1e-7,(filename,"T")
        assert np.all(a["R"]>=0) and np.all(a["T"]>=0)
        assert np.all(a["R"]+a["T"]<=1+1e-7)
        if not absorbing:
            assert np.max(np.abs(a["R"]+a["T"]-1))<1e-7
