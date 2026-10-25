"""Generalization checks across film thickness, dispersion and measurement noise.

Synthetic n/k models are explicitly identified as illustrative, not literature
measurements. Malitson (1965) provides the silica substrate dispersion.
"""
import numpy as np
import pytest
from engine import FitSettings,fit_absorption_aware_blind
from optical_inversion_package.inversion import Spectrum,coherent_film

def silica(w):
    l=np.asarray(w)/1000.
    return np.sqrt(1+.6961663*l*l/(l*l-.0684043**2)+.4079426*l*l/(l*l-.1162414**2)+.8974794*l*l/(l*l-9.896161**2))

@pytest.mark.parametrize('thickness_nm,dispersion,noise',[
    (900.,'cauchy',0.),
    (2000.,'cauchy',.002),
    (3000.,'cauchy',0.),
    (2000.,'absorbing',0.),
])
def test_blind_generalization(thickness_nm,dispersion,noise):
    w=np.linspace(450.,1200.,76)
    n=1.72+.085*(800./w)**2
    k=.012+.045*np.exp(-((w-510.)/105.)**2)
    if dispersion=='absorbing':
        n=1.95+.045*(800./w)**2
        k=.025+.09*np.exp(-((w-540.)/130.)**2)
    R,T=coherent_film(w,n+1j*k,thickness_nm,silica(w))
    if noise:
        rng=np.random.default_rng(20261025)
        R=np.clip(R+rng.normal(0,noise,len(w)),0,1)
        T=np.clip(T+rng.normal(0,noise,len(w)),0,1)
    spec=Spectrum(w,R,T,thickness_nm,substrate_n=silica(w))
    fit=fit_absorption_aware_blind(spec,FitSettings(n_basis=10,smoothness=.01,max_nfev=120),transparent_min_nm=750.,starts=5)
    n_mae=float(np.mean(np.abs(fit['n']-n)))
    k_mae=float(np.mean(np.abs(fit['k']-k)))
    assert fit['R_rmse']<.035, f'R RMSE={fit["R_rmse"]}, n MAE={n_mae}'
    assert fit['T_rmse']<.035, f'T RMSE={fit["T_rmse"]}, k MAE={k_mae}'
    assert n_mae<.12, f'n MAE={n_mae}'
    assert k_mae<.045, f'k MAE={k_mae}'
