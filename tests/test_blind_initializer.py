"""Regression checks for blind initialization on known synthetic optical constants."""
import numpy as np
import pytest
from engine import FitSettings,fit_absorption_aware_blind
from optical_inversion_package.inversion import Spectrum,coherent_film

def fused_silica(w):
    l=np.asarray(w)/1000.
    return np.sqrt(1+.6961663*l*l/(l*l-.0684043**2)+.4079426*l*l/(l*l-.1162414**2)+.8974794*l*l/(l*l-9.896161**2))

def test_blind_synthetic_known_nk():
    w=np.linspace(450.,1200.,76)
    n=1.72+.085*(800./w)**2
    k=.012+.045*np.exp(-((w-510.)/105.)**2)
    d=2000.
    R,T=coherent_film(w,n+1j*k,d,fused_silica(w))
    spec=Spectrum(w,R,T,d,substrate_n=fused_silica(w))
    result=fit_absorption_aware_blind(spec,FitSettings(n_basis=10,smoothness=.01,max_nfev=100),transparent_min_nm=750.,starts=2)
    assert np.all(np.isfinite(result['n'])) and np.all(np.isfinite(result['k']))
    assert np.all(result['k']>=-1e-8)
    assert result['R_rmse']<.035 and result['T_rmse']<.035
    assert np.mean(np.abs(result['n']-n))<.12
    assert np.mean(np.abs(result['k']-k))<.045

def test_blind_requires_both_channels():
    w=np.linspace(500.,1000.,25)
    spec=Spectrum(w,np.full(w.shape,.1),None,2000.,substrate_n=1.46)
    with pytest.raises(ValueError,match='both R and T'):
        fit_absorption_aware_blind(spec)
