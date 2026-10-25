import json
import numpy as np
from pathlib import Path
from optical_inversion_package.inversion import Spectrum,coherent_film
from engine import estimate_broadband_seed,phase2_hybrid_seed,fit_historical_logk
from historical_reference import archived_phase2

ROOT=Path(__file__).resolve().parent

def run():
    raw=json.loads((ROOT/'data/original_plot_raw_7.json').read_text())
    series={x['label']:x for x in raw['data']['series']}
    w=np.asarray(series['T']['x'],float)*1000
    R=np.asarray(series['R']['y'],float)
    T=np.asarray(series['T']['y'],float)
    from phase1_fresnel import phase1_fit
    mask=(w>=690)&(w<=1000)
    phase1=phase1_fit(w[mask],R[mask],T[mask],starts_nm=[1000],max_nfev=300,reference_nm=800)
    s=Spectrum(wavelength_nm=w,R=R,T=T,thickness_nm=phase1['thickness_nm'])
    x=800/w;p=phase1['parameters']
    n0=p['a']+p['b']*x*x+p['c']*x**4
    k0=np.exp(np.clip(p['q0']+p['q1']*(x-1),-24,2))
    broad=estimate_broadband_seed(s,n0,k0,max_nfev=70)
    hybrid=phase2_hybrid_seed(w,w,broad['n'],broad['k'],p)
    fit=fit_historical_logk(s,max_nfev=90,warm_start=hybrid)
    archive=archived_phase2()
    result=dict(phase1_thickness_nm=phase1['thickness_nm'],phase1_rmse=phase1['rmse'],broadband_nfev=broad['nfev'],phase2_nfev=fit['nfev'],phase2_candidates=fit['candidate_count'],phase2_success=fit['success'],phase2_train_rmse=fit['train_rmse'],phase2_holdout_rmse=fit['holdout_rmse'],archived_train_rmse=archive['train_rmse'],archived_holdout_rmse=archive['holdout_rmse'],n_mae=float(np.mean(np.abs(fit['n']-archive['n']))),k_mae=float(np.mean(np.abs(fit['k']-archive['k']))),thickness_difference_nm=phase1['thickness_nm']-archive['thickness_nm'])
    return result

if __name__=='__main__':
    print(json.dumps(run(),indent=2))
