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
    # The oracle-seeded test below measures optimizer accuracy only.
    # Independently test initialization without supplying true n/k.
    from engine import estimate_broadband_seed,phase2_hybrid_seed
    generic_n=np.full_like(w,1.7)
    generic_k=np.full_like(w,.03)
    broad=estimate_broadband_seed(s,generic_n,generic_k,nodes=25,max_nfev=70)
    p=dict(a=1.7,b=0.,c=0.,q0=float(np.log(.03)),q1=0.)
    hybrid=phase2_hybrid_seed(w,w,broad['n'],broad['k'],p,nodes=25)
    blind=fit_historical_logk(s,nodes=25,max_nfev=85,warm_start=hybrid)
    fit=fit_historical_logk(s,nodes=25,max_nfev=85,n_seed=n,k_seed=k)
    return dict(dataset='synthetic absorbing dye-like thin film (not literature experiment)',
        thickness_nm=d,noise_sigma=.001,phase2_success=fit['success'],
        training_rmse=fit['train_rmse'],holdout_rmse=fit['holdout_rmse'],
        n_mae=float(np.mean(np.abs(fit['n']-n))),
        k_mae=float(np.mean(np.abs(fit['k']-k))),
        candidate_count=fit['candidate_count'],nfev=fit['nfev'],\n        blind_training_rmse=blind['train_rmse'],blind_holdout_rmse=blind['holdout_rmse'],\n        blind_n_mae=float(np.mean(np.abs(blind['n']-n))),blind_k_mae=float(np.mean(np.abs(blind['k']-k))),\n        blind_success=blind['success'],blind_candidate_count=blind['candidate_count'])

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
