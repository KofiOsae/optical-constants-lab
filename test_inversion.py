import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from inversion import coherent_film, Spectrum, fit_joint, joint_information

def test_energy_conservation():
    wl=np.linspace(350,900,61)
    for angle in (0,30,60):
        for pol in ('s','p'):
            R,T=coherent_film(wl,np.full(len(wl),1.8+0j),200,1.5,angle_deg=angle,polarization=pol)
            np.testing.assert_allclose(R+T,1,atol=1e-12)

def test_zero_contrast():
    wl=np.linspace(400,800,30)
    R,T=coherent_film(wl,np.ones(len(wl),complex),200)
    np.testing.assert_allclose(R,0,atol=1e-12)
    np.testing.assert_allclose(T,1,atol=1e-12)

def test_passivity():
    wl=np.linspace(400,800,30)
    R,T=coherent_film(wl,np.full(len(wl),1.8+.1j),100,1.5)
    assert np.all(R>=0) and np.all(T>=0) and np.all(R+T<=1+1e-12)

def test_joint_recovery_synthetic():
    wl=np.linspace(400,800,65)
    truth=1.8+.025j
    spectra=[]
    for d,ns in [(60,1),(1000,1),(60,1.5),(1000,1.5)]:
        R,T=coherent_film(wl,np.full(len(wl),truth),d,ns)
        spectra.append(Spectrum(wl,R,T,d,ns))
    fit,fn=fit_joint(spectra,n_basis=8,n_init=1.75,k_init=.02,max_nfev=100)
    n,k=fn(wl)
    assert np.max(np.abs(n-truth.real))<.005
    assert np.max(np.abs(k-truth.imag))<.005

def test_information_improves_with_additional_measurements():
    a=[dict(thickness_nm=1000)]
    b=a+[dict(thickness_nm=50)]
    sa,_=joint_information(a,500,1.5,.7)
    sb,_=joint_information(b,500,1.5,.7)
    assert sb[-1]>=sa[-1]
