"""Dispersive multilayer joint inversion, isolated from existing single-film API.

Known layers are wavelength-aligned with each StackSpectrum. One unknown
material has shared dispersion; its thickness can vary by sample.
"""
import numpy as np
from scipy.optimize import least_squares
from optical_inversion_package.inversion import spline_design
from scientific_diagnostics import dielectric, jacobian_report
from multilayer_inverse import predict_stack

def fit_dispersive_stacks(samples, model='spline', n_basis=8, smoothness=.01,
                          oscillators=1, fit_thickness=False,
                          thickness_bounds=None, max_nfev=250, optical_bounds=None):
    if not samples: raise ValueError('At least one sample is required')
    if model not in ('spline','lorentz','drude_lorentz','kk_constrained'):
        raise ValueError('Unsupported dispersion model')
    if model!='spline' and not 1<=oscillators<=8:
        raise ValueError('Choose 1 to 3 oscillators')
    wl=np.concatenate([np.asarray(s.wavelength_nm,float) for s,d in samples])
    if np.any(~np.isfinite(wl)) or np.any(wl<=0):raise ValueError('Invalid wavelengths')
    lowwl,highwl=float(wl.min()),float(wl.max())
    if highwl<=lowwl:raise ValueError('Spectral range must be nonzero')
    nom=np.array([float(d) for s,d in samples])
    if np.any(nom<=0):raise ValueError('Positive thicknesses required')
    if optical_bounds is None: optical_bounds=(.5,6.,0.,5.)
    nlo,nhi,klo,khi=map(float,optical_bounds)
    if not np.all(np.isfinite([nlo,nhi,klo,khi])) or not (0<nlo<nhi and 0<=klo<khi):
        raise ValueError('Optical bounds require 0 < n_min < n_max and 0 <= k_min < k_max')
    if model=='spline':
        if n_basis<4:raise ValueError('Cubic spline needs at least 4 basis functions')
        B=[spline_design(s.wavelength_nm,lowwl,highwl,n_basis) for s,d in samples]
        x=np.r_[np.full(n_basis,1.8),np.full(n_basis,.02)]
        lo=np.r_[np.full(n_basis,nlo),np.full(n_basis,klo)]
        hi=np.r_[np.full(n_basis,nhi),np.full(n_basis,khi)]
        def nk(w,p,b=None):
            basis=b if b is not None else spline_design(w,lowwl,highwl,n_basis)
            return basis@p[:n_basis],basis@p[n_basis:2*n_basis]
    elif model=='kk_constrained':
        # Fixed positive Lorentz spectral basis: analytic, passive, KK-consistent.
        # Resonances span beyond measured band to mitigate truncation artifacts.
        elo,ehi=1239.841984/highwl,1239.841984/lowwl
        centers=np.geomspace(max(.02,elo*.35),ehi*2.5,oscillators)
        widths=np.maximum(.08*centers,.03)
        x=np.r_[2.5,np.full(oscillators,.25)]
        lo=np.r_[1.,np.zeros(oscillators)]
        hi=np.r_[25.,np.full(oscillators,30.)]
        B=[None]*len(samples)
        def nk(w,p,b=None):
            E=1239.841984/np.asarray(w,float)
            eps=np.full(E.shape,p[0],complex)
            for strength,center,width in zip(p[1:1+oscillators],centers,widths):
                eps+=strength*center**2/(center**2-E**2-1j*width*E)
            z=np.sqrt(eps+0j)
            return z.real,np.maximum(z.imag,0.)
    else:
        elo,ehi=1239.841984/highwl,1239.841984/lowwl
        x=[2.5];lo=[1.];hi=[25.]
        if model=='drude_lorentz':
            x.extend([.3,.15]);lo.extend([0.,.001]);hi.extend([10.,5.])
        for e in np.linspace(elo,ehi,oscillators+2)[1:-1]:
            x.extend([.3,float(e),.2])
            lo.extend([0.,.05,.001]);hi.extend([20.,max(10.,2*ehi),5.])
        x=np.asarray(x);lo=np.asarray(lo);hi=np.asarray(hi)
        B=[None]*len(samples)
        def nk(w,p,b=None):
            z=np.sqrt(dielectric(w,p[:disp_count],model,oscillators)+0j)
            return z.real,np.maximum(z.imag,0.)
    disp_count=len(x)
    if fit_thickness:
        if thickness_bounds is None or len(thickness_bounds)!=len(samples):
            raise ValueError('Provide per-sample thickness bounds')
        lower=np.array([v[0] for v in thickness_bounds],float)
        upper=np.array([v[1] for v in thickness_bounds],float)
        if np.any(lower<=0) or np.any(lower>=upper):
            raise ValueError('Invalid thickness bounds')
        x=np.r_[x,nom];lo=np.r_[lo,lower];hi=np.r_[hi,upper]
    def residual(p):
        chunks=[]
        for i,((s,d),b) in enumerate(zip(samples,B)):
            n,k=nk(np.asarray(s.wavelength_nm,float),p,b)
            thickness=p[disp_count+i] if fit_thickness else d
            r,t=predict_stack(s,n,k,thickness)
            for measured,predicted,sigma in ((s.R,r,s.sigma_R),(s.T,t,s.sigma_T)):
                if measured is not None:
                    obs=np.asarray(measured,float)
                    if obs.shape!=r.shape or not np.all(np.isfinite(obs)) or sigma<=0:
                        raise ValueError('Invalid measurements or uncertainties')
                    chunks.append((predicted-obs)/sigma)
        if model=='spline' and smoothness>0:
            chunks.extend([np.sqrt(smoothness)*np.diff(p[:n_basis],2),
                           np.sqrt(smoothness)*np.diff(p[n_basis:disp_count],2)])
        return np.concatenate(chunks)
    if model=='spline': x=np.clip(x,lo+1e-8,hi-1e-8)
    if not np.all((x>=lo)&(x<=hi)):raise ValueError('Initial parameters outside bounds')
    # Deterministic metal-aware multistart. The same objective and bounds are
    # used for every start; select by data+regularization residual norm.
    starts=[x]
    if model=='spline':
        for n_seed,k_seed in ((.5,2.),(1.,3.),(2.,2.5)):
            trial=x.copy()
            trial[:n_basis]=n_seed
            trial[n_basis:disp_count]=k_seed
            starts.append(np.clip(trial,lo+1e-8,hi-1e-8))
    elif model=='drude_lorentz':
        trial=x.copy()
        trial[1]=min(5.,hi[1]*.5)
        trial[2]=.15
        starts.append(trial)
    fits=[]
    for seed in starts:
        candidate=least_squares(residual,seed,bounds=(lo,hi),
                                max_nfev=max_nfev,x_scale='jac')
        fits.append(candidate)
    result=min(fits,key=lambda fit: float(np.dot(fit.fun,fit.fun)))
    def optical_constants(w):
        return nk(np.asarray(w,float),result.x)
    return {'result':result,'optical_constants':optical_constants,
            'thicknesses':result.x[disp_count:] if fit_thickness else nom,
            'weighted_rmse':float(np.sqrt(np.mean(result.fun**2))),
            'sensitivity':jacobian_report(result),
            'multistart':{'attempts':len(fits),'costs':[float(f.cost) for f in fits], 'selected_start':int(np.argmin([f.cost for f in fits]))},
            'optical_bounds':{'n_min':nlo,'n_max':nhi,'k_min':klo,'k_max':khi} if model=='spline' else None,
            'model':model,'dispersion_parameters':result.x[:disp_count].tolist(),
            'model_note':('Passive analytic fixed-pole Lorentz basis; KK consistency follows from the model, not from finite-band numerical integration. Basis is restrictive and must be tested for adequacy.' if model=='kk_constrained' else '')}
