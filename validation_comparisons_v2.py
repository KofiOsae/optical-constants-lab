"""Comparable, conservative inverse-solution diagnostics.

Inputs are actual repeated-fit records; never manufacture independent trials.
Outputs are descriptive screening metrics, not statistical confidence levels.
"""
import numpy as np

def compare_fits(runs, reference=None):
    if len(runs)<2: return {'status':'insufficient_runs','warning':'At least two completed fits required'}
    rows=[]
    for r in runs:
        p=np.asarray(r['parameters'],float)
        n=np.asarray(r['n'],float);k=np.asarray(r['k'],float)
        if not(len(n)==len(k)>0 and np.isfinite(p).all() and np.isfinite(n).all() and np.isfinite(k).all()):
            raise ValueError('Nonfinite or mismatched run arrays')
        rows.append((float(r['weighted_rmse']),p,n,k))
    if any(len(x[1])!=len(rows[0][1]) or len(x[2])!=len(rows[0][2]) for x in rows):
        raise ValueError('Compare aligned wavelengths and identical model families')
    best=min(rows,key=lambda x:x[0])
    rmse=np.array([x[0] for x in rows])
    spectra=np.array([np.r_[x[2],x[3]] for x in rows])
    spectral_spread=float(np.max(np.sqrt(np.mean((spectra-spectra.mean(axis=0))**2,axis=1))))
    out={'status':'evaluated','runs':len(rows),'best_weighted_rmse':float(rmse.min()),
         'rmse_max_over_min':float(rmse.max()/max(rmse.min(),1e-12)),
         'max_spectral_deviation_from_mean':spectral_spread,
         'parameter_ranges':(np.max([x[1] for x in rows],axis=0)-np.min([x[1] for x in rows],axis=0)).tolist(),
         'warning':'A stable local fit does not establish global uniqueness; same model, aligned wavelength grid and identical noise weights required.'}
    if reference is not None:
        nr=np.asarray(reference['n'],float);kr=np.asarray(reference['k'],float)
        if nr.shape!=best[2].shape or kr.shape!=best[3].shape: raise ValueError('Reference wavelength alignment required')
        out['reference_rmse_n']=float(np.sqrt(np.mean((best[2]-nr)**2)))
        out['reference_rmse_k']=float(np.sqrt(np.mean((best[3]-kr)**2)))
    return out

def window_stability(runs):
    """Repeated wavelength-window inversions, each with overlap-sampled fitted n,k."""
    out=compare_fits(runs)
    if out['status']=='evaluated':out['test']='window_stability'
    return out

def thickness_profile(runs):
    """Profile from actual thickness-perturbed fit records, not synthetic proxy."""
    out=compare_fits(runs)
    if out['status']=='evaluated':
        out['thicknesses_nm']=[float(x['thickness_nm']) for x in runs]
        out['test']='thickness_sensitivity'
    return out

def holdout_metrics(observed_R,observed_T,predicted_R,predicted_T):
    pairs=[(observed_R,predicted_R),(observed_T,predicted_T)]
    errors=[]
    for actual,pred in pairs:
        if actual is None or pred is None:continue
        a=np.asarray(actual,float);b=np.asarray(pred,float)
        if a.shape!=b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():raise ValueError('Invalid aligned holdout arrays')
        errors.extend((a-b).tolist())
    return {'available':bool(errors),'rmse':float(np.sqrt(np.mean(np.square(errors)))) if errors else None,
            'heldout_points':len(errors),
            'note':'Out-of-sample only if holdout was excluded from all fitting, initialization and tuning.'}
