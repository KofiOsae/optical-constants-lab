"""Synthetic physical recovery benchmark for the inversion engine.

Tests known n(lambda), k(lambda), thickness against generated R/T, with
a deterministic perturbation. This is not experimental literature validation.
"""
import json
import numpy as np
from optical_inversion_package.inversion import Spectrum,coherent_film
from engine import fit_historical_logk

def run():
    w=np.linspace(350.,950.,121)
    x=700/w
    n=1.62+.065*x*x
    k=.012+.18*np.exp(-.5*((w-530.)/65.)**2)
    d=820.
    R,T=coherent_film(w,n+1j*k,d)
    rng=np.random.default_rng(20261025)
    R_obs=np.clip(R+rng.normal(0,.001,len(w)),0,1)
    T_obs=np.clip(T+rng.normal(0,.001,len(w)),0,1)
    s=Spectrum(wavelength_nm=w,R=R_obs,T=T_obs,thickness_nm=d)
    fit=fit_historical_logk(s,nodes=25,max_nfev=85,n_seed=n,k_seed=k)
    return dict(dataset='synthetic absorbing dye-like thin film (not literature experiment)',
        thickness_nm=d,noise_sigma=.001,phase2_success=fit['success'],
        training_rmse=fit['train_rmse'],holdout_rmse=fit['holdout_rmse'],
        n_mae=float(np.mean(np.abs(fit['n']-n))),
        k_mae=float(np.mean(np.abs(fit['k']-k))),
        candidate_count=fit['candidate_count'],nfev=fit['nfev'])

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
