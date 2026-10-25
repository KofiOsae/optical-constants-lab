"""Thickness sensitivity for measured single-film R/T. No experimental data is stored here.
Usage: python verify_thickness_sensitivity.py RT_wavelength.csv --seed blind_initializer_result.csv
Seed is optional; with --seed this is a warm-start sensitivity diagnostic, NOT blind validation.
"""
import argparse
import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from optical_inversion_package.inversion import coherent_film,spline_design

def silica(w):
    l=np.asarray(w)/1000.
    return np.sqrt(1+.6961663*l*l/(l*l-.0684043**2)+.4079426*l*l/(l*l-.1162414**2)+.8974794*l*l/(l*l-9.896161**2))

def run(csv,seed=None,thicknesses=(1700,1800,1900,1950,2000,2050,2100,2200,2300),max_nfev=70):
    raw=pd.read_csv(csv)
    w=raw['wavelength_um'].to_numpy()*1000 if 'wavelength_um' in raw else raw['wavelength_nm'].to_numpy()
    R=raw['R'].to_numpy();T=raw['T'].to_numpy()
    m=18;B=spline_design(w,w[0],w[-1],m)
    lo=np.r_[np.ones(m),np.zeros(m)];hi=np.r_[np.full(m,4.),np.full(m,2.)]
    x0=np.r_[np.full(m,1.75),np.full(m,.05)]
    if seed is not None:
        prior=pd.read_csv(seed)
        sw=prior['wavelength_nm'].to_numpy()
        x0=np.r_[np.linalg.lstsq(B,np.interp(w,sw,prior['n']),rcond=None)[0],
                 np.linalg.lstsq(B,np.interp(w,sw,prior['k']),rcond=None)[0]]
    x0=np.clip(x0,lo+1e-8,hi-1e-8)
    rows=[]
    for d in thicknesses:
        def residual(x):
            rr,tt=coherent_film(w,B@x[:m]+1j*(B@x[m:]),d,silica(w))
            return np.r_[(rr-R)/.005,(tt-T)/.005,.1*np.diff(x[:m],2),.1*np.diff(x[m:],2)]
        fit=least_squares(residual,x0,bounds=(lo,hi),max_nfev=max_nfev)
        rr,tt=coherent_film(w,B@fit.x[:m]+1j*(B@fit.x[m:]),d,silica(w))
        rows.append({'thickness_nm':d,'R_rmse':float(np.sqrt(np.mean((rr-R)**2))),
            'T_rmse':float(np.sqrt(np.mean((tt-T)**2))),
            'combined_rmse':float(np.sqrt(np.mean(np.r_[(rr-R),(tt-T)]**2))),
            'n_800':float(np.interp(800,w,B@fit.x[:m])),
            'k_800':float(np.interp(800,w,B@fit.x[m:])),
            'nfev':fit.nfev,'success':bool(fit.success)})
    return pd.DataFrame(rows)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('csv');p.add_argument('--seed');p.add_argument('--output',default='thickness_sensitivity_results.csv');a=p.parse_args()
    out=run(a.csv,a.seed);out.to_csv(a.output,index=False);print(out.to_string(index=False))
