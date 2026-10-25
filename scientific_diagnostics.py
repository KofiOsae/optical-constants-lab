"""Causal Lorentz / Drude-Lorentz dielectric fits for coherent-film R/T.
Conventions: exp(-i omega t), N=n+i k, energy E=hc/lambda in eV.
"""
import numpy as np
from scipy.optimize import least_squares
from engine import coherent_film

def dielectric(wavelength_nm, params, model='lorentz', oscillators=1):
    E=1239.841984/np.asarray(wavelength_nm,float)
    eps=np.full(E.shape,params[0],dtype=complex)
    j=1
    if model=='drude_lorentz':
        plasma,gamma=params[1:3]; eps-=plasma**2/(E**2+1j*gamma*E);j=3
    for _ in range(oscillators):
        strength,center,width=params[j:j+3];j+=3
        eps+=strength*center**2/(center**2-E**2-1j*width*E)
    return eps

def fit_causal(spectra,model='lorentz',oscillators=1,max_nfev=180,fit_thickness=False):
    if model not in ('lorentz','drude_lorentz') or not 1<=oscillators<=3:
        raise ValueError('Choose 1–3 Lorentz oscillators, optionally with Drude term')
    allw=np.concatenate([s.wavelength_nm for s in spectra]);elo=1239.841984/max(allw);ehi=1239.841984/min(allw)
    x=[2.5];lo=[1.0];hi=[25.]
    if model=='drude_lorentz':x.extend([.3,.15]);lo.extend([0.,.001]);hi.extend([10.,5.])
    for e in np.linspace(elo,ehi,oscillators+2)[1:-1]:
        x.extend([.3,float(e),.2]);lo.extend([0.,.05,.001]);hi.extend([20.,max(10.,2*ehi),5.])
    nom=np.array([s.thickness_nm for s in spectra])
    if fit_thickness:
        x.extend(nom);lo.extend(nom*.95);hi.extend(nom*1.05)
    x=np.array(x);lo=np.array(lo);hi=np.array(hi)
    def optical(w,p):
        N=np.sqrt(dielectric(w,p,model,oscillators)+0j)
        return np.real(N),np.maximum(np.imag(N),0.)
    def residual(p):
        out=[]
        for i,s in enumerate(spectra):
            n,k=optical(s.wavelength_nm,p);d=p[-len(spectra)+i] if fit_thickness else s.thickness_nm
            R,T=coherent_film(s.wavelength_nm,n+1j*k,d,s.substrate_n,s.ambient_n,s.angle_deg,s.polarization)
            if s.R is not None:out.extend((R-s.R)/s.sigma_R)
            if s.T is not None:out.extend((T-s.T)/s.sigma_T)
        return np.array(out)
    result=least_squares(residual,x,bounds=(lo,hi),max_nfev=max_nfev,x_scale='jac')
    p=result.x.copy()
    return dict(result=result,optical_constants=lambda w:optical(w,p),
        thicknesses=p[-len(spectra):] if fit_thickness else nom,
        weighted_rmse=float(np.sqrt(np.mean(residual(p)**2))),
        stages={'selected':model,'oscillators':oscillators,'coefficients':p.tolist()},
        model=model)

def jacobian_report(result, noise_is_known=True):
    J=np.asarray(result.jac,float)
    if J.ndim!=2 or not np.all(np.isfinite(J)):return {'available':False}
    u,s,v=np.linalg.svd(J,full_matrices=False)
    threshold=max(J.shape)*np.finfo(float).eps*s[0] if len(s) else 0
    rank=int(np.sum(s>threshold))
    # Covariance is local linear approximation; do not interpret as global uniqueness.
    cov=np.linalg.pinv(J.T@J,rcond=1e-12)
    std=np.sqrt(np.maximum(0,np.diag(cov)))
    return {'available':True,'rank':rank,'parameters':J.shape[1],
      'observations':J.shape[0],'condition_number':float(s[0]/s[-1]) if len(s) and s[-1]>0 else None,
      'singular_values':s.tolist(),'parameter_standard_errors_local':std.tolist(),
      'note':'Local linearized sensitivity; assumes independent, correctly normalized measurement errors. Does not prove global uniqueness. Regularization may affect the Jacobian.'}
