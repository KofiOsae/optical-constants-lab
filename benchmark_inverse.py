import json,time
import numpy as np,pandas as pd
from scipy.optimize import least_squares,minimize_scalar
from engine import FitSettings,fit
from optical_inversion_package.inversion import Spectrum,coherent_film
from pathlib import Path
root=Path(__file__).resolve().parent
raw=json.loads((root/'data/original_plot_raw_7.json').read_text())
s={v['label']:v for v in raw['data']['series']}
w=np.array(s['R']['x'])*1000;R=np.array(s['R']['y']);T=np.array(s['T']['y']);mask=w>=690
spec=Spectrum(w[mask],R[mask],T[mask],1000.,name='Original R6G 1um')
start=time.time()
settings=FitSettings(n_basis=10,smoothness=.1,optimizer='least_squares',fit_thickness=True,thickness_tolerance=.15,max_nfev=180)
r=fit([spec],settings); n,k=r['optical_constants'](w[mask]);d=float(r['thicknesses'][0]);Rp,Tp=coherent_film(w[mask],n+1j*k,d)
print('blind spline d',d,'spectral RMSE',np.sqrt(np.mean(np.r_[Rp-R[mask],Tp-T[mask]]**2)),'elapsed',time.time()-start)
ref=pd.read_csv(root/'data/R6G_Singh.csv');wr=ref.wavelength_um.to_numpy()*1000
nr=np.interp(w[mask],wr,ref.n);kr=np.interp(w[mask],wr,ref.k)
def residual_d(d):
 a,b=coherent_film(w[mask],nr+1j*kr,d);return np.r_[a-R[mask],b-T[mask]]
sol=least_squares(lambda x:residual_d(x[0]),[980.],bounds=([850.],[1150.]),xtol=1e-13,ftol=1e-13,gtol=1e-13)
print('reference assisted d',sol.x[0],np.sqrt(np.mean(sol.fun**2)))
np.savez(root/'data/original_inverse_live_results.npz',wavelength_nm=w[mask],R_observed=R[mask],T_observed=T[mask],R_fit=Rp,T_fit=Tp,n_fit=n,k_fit=k,n_reference=nr,k_reference=kr,thickness_fit_nm=d,reference_assisted_thickness_nm=sol.x[0])
