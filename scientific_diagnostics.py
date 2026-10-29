"""Causal Lorentz / Drude-Lorentz dielectric fits for coherent-film R/T.
Conventions: exp(-i omega t), N=n+i k, energy E=hc/lambda in eV.
"""
import numpy as np
from scipy.optimize import least_squares
from engine import coherent_film
from organic_dispersion import MODEL_PARAMETER_COUNTS, organic_parameter_setup, organic_dielectric

def dielectric(wavelength_nm, params, model='lorentz', oscillators=1):
    if model in MODEL_PARAMETER_COUNTS:
        return organic_dielectric(wavelength_nm,params,model,oscillators)
    E=1239.841984/np.asarray(wavelength_nm,float)
    eps=np.full(E.shape,params[0],dtype=complex)
    j=1
    if model=='drude_lorentz':
        plasma,gamma=params[1:3]; eps-=plasma**2/(E**2+1j*gamma*E);j=3
    for _ in range(oscillators):
        strength,center,width=params[j:j+3];j+=3
        eps+=strength*center**2/(center**2-E**2-1j*width*E)
    return eps

def fit_causal(spectra,model='lorentz',oscillators=1,max_nfev=180,fit_thickness=False,multistart=1,random_seed=20261026,warm_start=None,parameter_controls=None):
    if model not in ('lorentz','drude_lorentz',*MODEL_PARAMETER_COUNTS) or not 1<=oscillators<=8:
        raise ValueError('Choose 1–3 Lorentz oscillators, optionally with Drude term')
    allw=np.concatenate([s.wavelength_nm for s in spectra]);elo=1239.841984/max(allw);ehi=1239.841984/min(allw)
    if model in MODEL_PARAMETER_COUNTS:
        x,lo,hi=organic_parameter_setup(model,oscillators,min(allw),max(allw))
        x=x.tolist();lo=lo.tolist();hi=hi.tolist()
    else:
        x=[2.5];lo=[1.0];hi=[25.]
    if model=='drude_lorentz':x.extend([.3,.15]);lo.extend([0.,.001]);hi.extend([10.,5.])
    if model not in MODEL_PARAMETER_COUNTS:
        for e in np.linspace(elo,ehi,oscillators+2)[1:-1]:
            x.extend([.3,float(e),.2]);lo.extend([0.,.05,.001]);hi.extend([20.,max(10.,2*ehi),5.])
    nom=np.array([s.thickness_nm for s in spectra])
    if fit_thickness:
        x.extend(nom);lo.extend(nom*.95);hi.extend(nom*1.05)
    x=np.array(x);lo=np.array(lo);hi=np.array(hi)
    if parameter_controls is not None:
        if not isinstance(parameter_controls,dict):
            raise ValueError('Parameter controls must be an object')
        count=len(x)-(len(spectra) if fit_thickness else 0)
        for key,arr in [('initial',x),('lower',lo),('upper',hi)]:
            if key in parameter_controls:
                values=np.asarray(parameter_controls[key],float)
                if values.shape!=(count,) or not np.all(np.isfinite(values)):
                    raise ValueError(f'{key} requires {count} finite optical coefficients')
                arr[:count]=values
        if np.any(lo>=hi) or np.any(x<lo) or np.any(x>hi):
            raise ValueError('Initial parameters must lie within strict bounds')

    def optical(w,p):
        N=np.sqrt(dielectric(w,p,model,oscillators)+0j)
        return np.real(N),np.maximum(np.imag(N),0.)
    def residual(p):
        out=[]
        for i,s in enumerate(spectra):
            n,k=optical(s.wavelength_nm,p);d=p[-len(spectra)+i] if fit_thickness else s.thickness_nm
            R,T=coherent_film(s.wavelength_nm,n+1j*k,d,s.substrate_n,s.ambient_n,s.angle_deg,s.polarization,roughness_nm=getattr(s,'roughness_nm',0.))
            if s.R is not None:out.extend((R-s.R)/s.sigma_R)
            if s.T is not None:out.extend((T-s.T)/s.sigma_T)
        return np.array(out)
    if not isinstance(multistart,int) or not 1 <= multistart <= 50:
        raise ValueError('multistart must be an integer from 1 to 50')
    if warm_start is not None:
        previous=np.asarray(warm_start,dtype=float)
        offset=3 if model=='drude_lorentz' else 1
        stride=MODEL_PARAMETER_COUNTS.get(model,3)
        expected=offset+stride*(oscillators-1)
        if fit_thickness: raise ValueError('Warm start requires fixed film thicknesses')
        if previous.shape!=(expected,) or not np.all(np.isfinite(previous)):
            raise ValueError('Warm-start coefficients must represent the same causal model with one fewer oscillator')
        x[:expected]=previous
        x[expected:expected+stride]=[1e-6,*x[expected+1:expected+stride]]
    rng=np.random.default_rng(random_seed)
    candidates=[]
    for run in range(multistart):
        seed=x.copy()
        if run:
            seed[0]=rng.uniform(1.0,6.0)
            offset=3 if model=='drude_lorentz' else 1
            if model in MODEL_PARAMETER_COUNTS:
                for j in range(oscillators):
                    idx=offset+MODEL_PARAMETER_COUNTS[model]*j
                    for k in range(MODEL_PARAMETER_COUNTS[model]):
                        seed[idx+k]=rng.uniform(lo[idx+k]+.01*(hi[idx+k]-lo[idx+k]),
                                                  lo[idx+k]+.45*(hi[idx+k]-lo[idx+k]))
            if model=='drude_lorentz':
                seed[1:3]=[rng.uniform(.05,2.0),rng.uniform(.02,1.0)]
            for j in range(oscillators if model not in MODEL_PARAMETER_COUNTS else 0):
                idx=offset+3*j
                seed[idx:idx+3]=[rng.uniform(.05,2.0),
                    rng.uniform(max(.05,elo*.75),min(hi[idx+1],ehi*1.25)),
                    rng.uniform(.02,.6)]
            if fit_thickness:
                seed[-len(spectra):]=nom*rng.uniform(.98,1.02,size=len(spectra))
        seed=np.clip(seed,lo+1e-9*(hi-lo),hi-1e-9*(hi-lo))
        starting_rmse=float(np.sqrt(np.mean(residual(seed)**2)))
        trial=least_squares(residual,seed,bounds=(lo,hi),max_nfev=max_nfev,x_scale='jac')
        candidates.append((float(np.mean(trial.fun**2)),trial,run+1,seed.tolist(),starting_rmse))
    candidates.sort(key=lambda item:item[0])
    result=candidates[0][1]
    p=result.x.copy()
    return dict(result=result,optical_constants=lambda w:optical(w,p),
        thicknesses=p[-len(spectra):] if fit_thickness else nom,
        weighted_rmse=float(np.sqrt(np.mean(residual(p)**2))),
        stages={'selected':model,'oscillators':oscillators,'coefficients':p.tolist(),'multistart':multistart,'random_seed':random_seed,'warm_start_used':warm_start is not None,'start_results':[{'initialization':int(start),'weighted_rmse':float(np.sqrt(score)),'success':bool(trial.success),'status':int(trial.status),'nfev':int(trial.nfev),'message':str(trial.message),'selected':bool(i==0),'coefficients':trial.x.tolist(),'initial_coefficients':seed,'initial_weighted_rmse':start_rmse} for i,(score,trial,start,seed,start_rmse) in enumerate(candidates)]},
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
    # Rank is a numerical property of the local Jacobian, not proof of unique n,k.
    # A conservative warning avoids universal condition-number pass/fail thresholds.
    smallest=float(s[-1]) if len(s) else 0.
    largest=float(s[0]) if len(s) else 0.
    relative_smallest=smallest/largest if largest>0 else 0.
    warnings=[]
    if rank<J.shape[1]: warnings.append('Local Jacobian is rank deficient: at least one parameter combination is not resolved.')
    if relative_smallest<1e-6: warnings.append('Very weak local sensitivity direction (smallest/largest singular value < 1e-6); inspect parameter stability.')
    if J.shape[0]<=J.shape[1]: warnings.append('Number of observations does not exceed the number of fitted parameters.')
    if not warnings: warnings.append('No severe numerical warning detected; global uniqueness and model adequacy remain unverified.')
    return {'available':True,'relative_smallest_singular_value':relative_smallest,
      'identifiability_warnings':warnings,'rank':rank,'parameters':J.shape[1],
      'observations':J.shape[0],'condition_number':float(s[0]/s[-1]) if len(s) and s[-1]>0 else None,
      'singular_values':s.tolist(),'parameter_standard_errors_local':std.tolist(),
      'note':'Local linearized sensitivity; assumes independent, correctly normalized measurement errors. Does not prove global uniqueness. Regularization may affect the Jacobian.'}
