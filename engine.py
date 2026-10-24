from dataclasses import dataclass
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import least_squares, differential_evolution
from optical_inversion_package.inversion import Spectrum,coherent_film,spline_design,local_sensitivity,joint_information

@dataclass
class FitSettings:
    n_basis:int=10
    smoothness:float=0.1
    optimizer:str='least_squares'
    fit_thickness:bool=False
    thickness_tolerance:float=.15
    max_nfev:int=180
    seed:int=42
    population:int=18
    iterations:int=30
    noise_floor:float=.005
    initial_n:object=None
    initial_thickness_nm:float=None

def validate_spectra(spectra):
    if not spectra: raise ValueError('At least one spectrum is required')
    for s in spectra:
        w=np.asarray(s.wavelength_nm)
        if w.ndim!=1 or len(w)<5 or np.any(np.diff(w)<=0): raise ValueError(f'{s.name}: wavelengths must be strictly increasing')
        if s.R is None and s.T is None: raise ValueError(f'{s.name}: missing R and T')
        if s.thickness_nm<=0: raise ValueError('Thickness must be positive')
        for kind in ('R','T'):
            arr=getattr(s,kind)
            if arr is not None and (len(arr)!=len(w) or not np.all(np.isfinite(arr)) or np.any((np.asarray(arr)<0)|(np.asarray(arr)>1))): raise ValueError(f'{s.name}: invalid {kind}')
        if s.sigma_R<=0 or s.sigma_T<=0: raise ValueError('Positive measurement uncertainties required')

def make_problem(spectra,settings):
    validate_spectra(spectra)
    m=settings.n_basis
    if m<4: raise ValueError('Need >=4 spline basis functions')
    lo=min(float(s.wavelength_nm[0]) for s in spectra); hi=max(float(s.wavelength_nm[-1]) for s in spectra)
    mats=[spline_design(s.wavelength_nm,lo,hi,m) for s in spectra]
    thicknesses=np.array([s.thickness_nm for s in spectra]); n_thick=len(spectra) if settings.fit_thickness else 0
    lower=np.r_[np.full(m,1.),np.zeros(m),thicknesses*(1-settings.thickness_tolerance) if n_thick else []]
    upper=np.r_[np.full(m,4.),np.full(m,2.),thicknesses*(1+settings.thickness_tolerance) if n_thick else []]
    x0=np.r_[np.full(m,1.75),np.full(m,.05),thicknesses if n_thick else []]
    if settings.initial_n is not None:
        seed_w, seed_n = settings.initial_n
        seed_w=np.asarray(seed_w,float); seed_n=np.asarray(seed_n,float)
        if len(seed_w)<2 or len(seed_w)!=len(seed_n) or np.any(np.diff(seed_w)<=0):
            raise ValueError('Invalid Swanepoel index seed')
        grid=np.linspace(lo,hi,max(4*m,50))
        B_seed=spline_design(grid,lo,hi,m)
        x0[:m]=np.clip(np.linalg.lstsq(B_seed,np.interp(grid,seed_w,seed_n),rcond=None)[0],low[:m]+1e-5,high[:m]-1e-5)
    if settings.initial_thickness_nm is not None and n_thick:
        d_seed=float(settings.initial_thickness_nm)
        if not np.isfinite(d_seed) or d_seed<=0:
            raise ValueError('Invalid Swanepoel thickness seed')
        x0[2*m:]=np.clip(d_seed,low[2*m:]+1e-5,high[2*m:]-1e-5)
    def residual(x,include_penalty=True):
        a=x[:m]; b=x[m:2*m]; d=x[2*m:] if n_thick else thicknesses
        out=[]
        for s,B,di in zip(spectra,mats,d):
            n=B@a; k=B@b
            R,T=coherent_film(s.wavelength_nm,n+1j*k,di,s.substrate_n,s.ambient_n,s.angle_deg,s.polarization)
            if s.R is not None: out.extend((R-s.R)/s.sigma_R)
            if s.T is not None: out.extend((T-s.T)/s.sigma_T)
        if include_penalty:
            out.extend(np.sqrt(settings.smoothness)*np.diff(a,2))
            out.extend(np.sqrt(settings.smoothness)*np.diff(b,2))
        return np.asarray(out)
    return residual,x0,lower,upper,(lo,hi),thicknesses

def harris_hawks(objective,lower,upper,population=18,iterations=30,seed=42):
    '''Bounded HHO: exploration, soft/hard besiege, and rapid dives with Levy steps.'''
    rng=np.random.default_rng(seed); dim=len(lower)
    hawks=rng.uniform(lower,upper,size=(population,dim))
    scores=np.array([objective(x) for x in hawks]); best=hawks[np.argmin(scores)].copy(); best_score=float(np.min(scores)); history=[best_score]; evaluations=population
    for t in range(iterations):
        energy0=2*rng.random()-1
        for i in range(population):
            E=2*energy0*(1-t/max(iterations,1)); q=rng.random(); r=rng.random()
            if abs(E)>=1:
                peer=hawks[rng.integers(population)]
                if q<.5: candidate=peer-rng.random(dim)*abs(peer-2*rng.random(dim)*hawks[i])
                else: candidate=(best-hawks.mean(axis=0))-rng.random(dim)*(lower+rng.random(dim)*(upper-lower))
            elif r>=.5:
                jump=2*(1-rng.random()); candidate=(best-E*abs(jump*best-hawks[i])) if abs(E)>=.5 else (best-E*abs(best-hawks[i]))
            else:
                jump=2*(1-rng.random()); center=hawks[i] if abs(E)>=.5 else hawks.mean(axis=0)
                candidate=best-E*abs(jump*best-center)
            candidate=np.clip(candidate,lower,upper); score=objective(candidate); evaluations+=1
            if score<scores[i]: hawks[i]=candidate; scores[i]=score
            elif abs(E)<1 and r<.5:
                # progressive rapid dive, heavy-tailed step
                levy=rng.standard_normal(dim)/np.maximum(abs(rng.standard_normal(dim)),1e-6)**(1/1.5)
                dive=np.clip(candidate+.01*(upper-lower)*levy,lower,upper)
                dive_score=objective(dive); evaluations+=1
                if dive_score<scores[i]: hawks[i]=dive; scores[i]=dive_score
            if scores[i]<best_score: best=hawks[i].copy(); best_score=float(scores[i])
        history.append(best_score)
    return best,np.asarray(history),evaluations

def fit(spectra,settings):
    residual,x0,low,high,domain,original_thicknesses=make_problem(spectra,settings)
    objective=lambda x: float(np.dot(residual(x),residual(x)))
    history=[]; evaluations_global=0
    if settings.optimizer=='least_squares': initial=x0
    elif settings.optimizer=='multistart':
        rng=np.random.default_rng(settings.seed); candidates=[x0]+[rng.uniform(low,high) for _ in range(4)]
        fits=[least_squares(residual,c,bounds=(low,high),max_nfev=max(25,settings.max_nfev//5)) for c in candidates]
        initial=min(fits,key=lambda z: np.linalg.norm(z.fun)).x
        history=[float(np.dot(f.fun,f.fun)) for f in fits]; evaluations_global=sum(f.nfev for f in fits)
    elif settings.optimizer=='differential_evolution':
        res=differential_evolution(objective,list(zip(low,high)),seed=settings.seed,popsize=max(3,settings.population//len(x0)),maxiter=settings.iterations,polish=False)
        initial=res.x; history=[float(res.fun)]; evaluations_global=res.nfev
    elif settings.optimizer in ('hho','hho_hybrid'):
        initial,history,evaluations_global=harris_hawks(objective,low,high,settings.population,settings.iterations,settings.seed)
    else: raise ValueError('Unknown optimizer')
    if settings.optimizer=='hho':
        class R: pass
        result=R(); result.x=initial; result.fun=residual(initial); result.nfev=evaluations_global; result.success=False; result.message='HHO search budget completed; numerical convergence not established (no local refinement)'
    else:
        result=least_squares(residual,initial,bounds=(low,high),max_nfev=settings.max_nfev,x_scale='jac')
    m=settings.n_basis
    def optical_constants(w):
        B=spline_design(np.asarray(w),*domain,m)
        return B@result.x[:m],B@result.x[m:2*m]
    thickness=result.x[2*m:] if settings.fit_thickness else original_thicknesses
    measurements=np.concatenate([residual(result.x,False)])
    return {'result':result,'optical_constants':optical_constants,'thicknesses':thickness,'weighted_rmse':float(np.sqrt(np.mean(measurements**2))),'history':history,'global_evaluations':evaluations_global,'domain':domain}

def simulate_four(wavelength_nm,n,k):
    specs=[]
    for name,d,sub in [('A thick air',1000.,1.),('B thin air',50.,1.),('C thick glass',1000.,1.5),('D thin glass',50.,1.5)]:
        R,T=coherent_film(wavelength_nm,n+1j*k,d,sub)
        specs.append(Spectrum(np.array(wavelength_nm),R,T,d,substrate_n=sub,name=name))
    return specs

def compare_prediction(spectrum,n,k,thickness=None):
    R,T=coherent_film(spectrum.wavelength_nm,np.asarray(n)+1j*np.asarray(k),spectrum.thickness_nm if thickness is None else thickness,spectrum.substrate_n,spectrum.ambient_n,spectrum.angle_deg,spectrum.polarization)
    diffs=np.concatenate([R-spectrum.R if spectrum.R is not None else np.array([]),T-spectrum.T if spectrum.T is not None else np.array([])])
    return {'R':R,'T':T,'rmse':float(np.sqrt(np.mean(diffs**2)))}
