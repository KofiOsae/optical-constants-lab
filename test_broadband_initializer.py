import numpy as np
from engine import estimate_broadband_seed,phase2_hybrid_seed
from optical_inversion_package.inversion import Spectrum,coherent_film

def test_sample_specific_broadband_hybrid_initialization():
    w=np.linspace(300.,1000.,41)
    n=np.full_like(w,1.72)
    k=np.full_like(w,.012)
    R,T=coherent_film(w,n+1j*k,1000.)
    s=Spectrum(name='synthetic',wavelength_nm=w,R=R,T=T,thickness_nm=1000.)
    broadband=estimate_broadband_seed(s,n,k,nodes=25,max_nfev=2)
    p=dict(a=1.72,b=0.,c=0.,q0=np.log(.012),q1=0.)
    q=phase2_hybrid_seed(w,w,broadband['n'],broadband['k'],p)
    assert q.shape==(50,)
    assert np.all(np.isfinite(q))
