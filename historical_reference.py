"""Immutable R6G regression fixture. Not a generic-material inversion seed."""
import json
from pathlib import Path
import numpy as np
from scipy.interpolate import CubicSpline
ROOT=Path(__file__).resolve().parent

def archived_phase2():
    fixture=json.loads((ROOT/'data/r6g_phase2_archived.json').read_text())
    raw=json.loads((ROOT/'data/original_plot_raw_7.json').read_text())
    series={item['label']:item for item in raw['data']['series']}
    w=np.asarray(series['T']['x'],dtype=float)
    T=np.asarray(series['T']['y'],dtype=float)
    R=np.asarray(series['R']['y'],dtype=float)
    # Workbench archive wavelength is in micrometres.
    q=np.asarray(fixture['spline'],dtype=float)
    knots=np.linspace(fixture['wavelength_min_nm'],fixture['wavelength_max_nm'],25)
    wn=w*1000
    n=CubicSpline(knots,q[:25],bc_type='natural')(wn)
    k=np.exp(np.clip(CubicSpline(knots,q[25:],bc_type='natural')(wn),-24,2))
    z=n+1j*k
    d=fixture['thickness_nm']
    r=(1-z)/(1+z)
    p=np.exp(2j*np.pi*z*d/wn)
    den=1-r*r*p*p
    tf=np.abs((1-r*r)*p/den)**2
    rf=np.abs(r*(1-p*p)/den)**2
    train=np.arange(len(w))%5!=0
    def rmse(mask):
        return float(np.sqrt(np.mean(np.r_[tf[mask]-T[mask],rf[mask]-R[mask]]**2)))
    return dict(wavelength_nm=wn.tolist(),n=n.tolist(),k=k.tolist(),R_fit=rf.tolist(),T_fit=tf.tolist(),thickness_nm=d,train_rmse=rmse(train),holdout_rmse=rmse(~train),coefficient_count=len(q),reference_only=True)
