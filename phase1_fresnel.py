"""Documented R6G Phase I: air/film/air coherent R/T inverse."""
import numpy as np
from scipy.optimize import least_squares
from optical_inversion_package.inversion import coherent_film

def phase1_fit(wavelength_nm, reflectance, transmittance, starts_nm=None, max_nfev=250, substrate_n=1.0, ambient_n=1.0, angle_deg=0.0, polarization='s'):
    w=np.asarray(wavelength_nm,float)
    R0=np.asarray(reflectance,float); T0=np.asarray(transmittance,float)
    if w.ndim!=1 or len(w)<8 or R0.shape!=w.shape or T0.shape!=w.shape:
        raise ValueError("Matching R/T wavelength arrays with at least 8 points required")
    if np.any(np.diff(w)<=0) or np.any(w<=0):
        raise ValueError('Phase I requires positive, strictly increasing wavelengths')
    if not np.isfinite(substrate_n) or substrate_n<1 or not np.isfinite(ambient_n) or ambient_n<=0:\n        raise ValueError('Phase I requires known positive real substrate and ambient indices')\n    x=float(np.median(w))/w
    def forward(p):
        a,b,c,q0,q1,d=p
        n=a+b*x*x+c*x**4
        k=np.exp(np.clip(q0+q1*(x-1),-35,3))
        R,T=coherent_film(w,n+1j*k,d,substrate_n,ambient_n,angle_deg,polarization)
        return n,k,R,T
    def residual(p):
        n,k,R,T=forward(p)
        return np.r_[R-R0,T-T0]
    lower=[1.,-3.,-3.,-30.,-40.,100.]
    upper=[4.,3.,3.,0.,40.,5000.]
    starts=starts_nm if starts_nm is not None else [450,600,750,900,1000,1100,1250,1500]
    trials=[]
    for d in starts:
        guess=[1.7375,-.0664,.1083,-9.23,12.56,float(d)]
        sol=least_squares(residual,guess,bounds=(lower,upper),max_nfev=max_nfev)
        rmse=float(np.sqrt(np.mean(residual(sol.x)**2)))
        trials.append((rmse,sol,float(d)))
    score,sol,start=min(trials,key=lambda item:item[0])
    n,k,R,T=forward(sol.x)
    return dict(geometry='freestanding' if abs(substrate_n-ambient_n)<1e-9 else 'film_on_semi_infinite_substrate',wavelength_nm=w.tolist(),n=n.tolist(),k=k.tolist(),
                R_fit=R.tolist(),T_fit=T.tolist(),R_measured=R0.tolist(),
                T_measured=T0.tolist(),thickness_nm=float(sol.x[5]),rmse=score,
                parameters=dict(zip(('a','b','c','q0','q1','d_nm'),map(float,sol.x))),
                success=bool(sol.success),message=str(sol.message),
                multistart=[dict(start_nm=d,thickness_nm=float(s.x[5]),rmse=v)
                            for v,s,d in trials])
