"""Reproduce original Phase II starting coefficients using archived global fit.

This is an R6G-only scientific regression. Never use the archived global
optical constants to initialize a different sample.
"""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from engine import phase2_hybrid_seed,fit_historical_logk
from historical_reference import archived_phase2
from optical_inversion_package.inversion import Spectrum

ROOT=Path(__file__).resolve().parent

def run():
    raw=json.loads((ROOT/'data/original_plot_raw_7.json').read_text())
    series={x['label']:x for x in raw['data']['series']}
    w=np.asarray(series['T']['x'],float)*1000
    R=np.asarray(series['R']['y'],float)
    T=np.asarray(series['T']['y'],float)
    global_fit=pd.read_csv(ROOT/'data/phase2_global_results.csv')
    p=dict(a=1.73753613,b=-.06640499,c=.10833076,q0=-9.230267,q1=12.56295895)
    warm=phase2_hybrid_seed(w,np.asarray(global_fit['wl_um'])*1000,
                           np.asarray(global_fit['n_global']),np.asarray(global_fit['k_global']),p)
    archive=archived_phase2()
    s=Spectrum(wavelength_nm=w,R=R,T=T,thickness_nm=archive['thickness_nm'])
    fit=fit_historical_logk(s,max_nfev=90,warm_start=warm)
    return dict(phase2_train_rmse=fit['train_rmse'],phase2_holdout_rmse=fit['holdout_rmse'],
                archived_train_rmse=archive['train_rmse'],archived_holdout_rmse=archive['holdout_rmse'],
                n_mae=float(np.mean(np.abs(fit['n']-archive['n']))),
                k_mae=float(np.mean(np.abs(fit['k']-archive['k']))),
                nfev=fit['nfev'],success=fit['success'])
if __name__=='__main__':
    print(json.dumps(run(),indent=2))
