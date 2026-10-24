import json
import numpy as np
import pandas as pd
from pathlib import Path
from engine import FitSettings,fit
from optical_inversion_package.inversion import Spectrum,coherent_film
from scipy.optimize import least_squares
ROOT=Path(__file__).resolve().parent

def run_original_inverse(initial_thickness_nm=950.0):
    if not 850 <= initial_thickness_nm <= 1150: raise ValueError('Initial thickness must be 850–1150 nm')
    raw=json.loads((ROOT/'data/original_plot_raw_7.json').read_text())
    curves={v['label']:v for v in raw['data']['series']}
    w=np.asarray(curves['R']['x'])*1000
    R=np.asarray(curves['R']['y']);T=np.asarray(curves['T']['y']);selection=w>=690
    w,R,T=w[selection],R[selection],T[selection]
    s=Spectrum(w,R,T,initial_thickness_nm,name='Original Workbench spectrum')
    result=fit([s],FitSettings(n_basis=10,smoothness=.1,fit_thickness=True,max_nfev=180))
    n,k=result['optical_constants'](w); d=float(result['thicknesses'][0]);rf,tf=coherent_film(w,n+1j*k,d)
    ref=pd.read_csv(ROOT/'data/R6G_Singh.csv')
    nref=np.interp(w,ref.wavelength_um.to_numpy()*1000,ref.n)
    kref=np.interp(w,ref.wavelength_um.to_numpy()*1000,ref.k)
    assisted=least_squares(lambda x:np.r_[coherent_film(w,nref+1j*kref,x[0])[0]-R,coherent_film(w,nref+1j*kref,x[0])[1]-T],[initial_thickness_nm],bounds=([800],[1200]),max_nfev=150)
    return {'wavelength_nm':w.tolist(),'R_observed':R.tolist(),'T_observed':T.tolist(),'R_fit':rf.tolist(),'T_fit':tf.tolist(),'n_fit':n.tolist(),'k_fit':k.tolist(),'n_reference':nref.tolist(),'k_reference':kref.tolist(),'initial_thickness_nm':initial_thickness_nm,'blind_spline_thickness_nm':d,'reference_assisted_thickness_nm':float(assisted.x[0]),'blind_spectral_rmse':float(np.sqrt(np.mean(np.r_[rf-R,tf-T]**2))),'reference_assisted_spectral_rmse':float(np.sqrt(np.mean(assisted.fun**2))),'n_rmse':float(np.sqrt(np.mean((n-nref)**2))),'k_rmse':float(np.sqrt(np.mean((k-kref)**2))),'converged':bool(result['result'].success),'note':'Blind spline fit uses a supplied initial thickness and bounds ±15%; a low spectral RMSE does not establish thickness uniqueness. Reference-assisted thickness uses known n,k and is not blind.'}
