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
    n_min:float=0.05
    n_max:float=8.0
    k_max:float=6.0

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
    lower=np.r_[np.full(m,settings.n_min),np.zeros(m),thicknesses*(1-settings.thickness_tolerance) if n_thick else []]
    upper=np.r_[np.full(m,settings.n_max),np.full(m,settings.k_max),thicknesses*(1+settings.thickness_tolerance) if n_thick else []]
    x0=np.r_[np.full(m,1.75),np.full(m,.05),thicknesses if n_thick else []]
    if settings.initial_n is not None:
        seed_w, seed_n = settings.initial_n
        seed_w=np.asarray(seed_w,float); seed_n=np.asarray(seed_n,float)
        if len(seed_w)<2 or len(seed_w)!=len(seed_n) or np.any(np.diff(seed_w)<=0):
            raise ValueError('Invalid Swanepoel index seed')
        grid=np.linspace(lo,hi,max(4*m,50))
        B_seed=spline_design(grid,lo,hi,m)
        x0[:m]=np.clip(np.linalg.lstsq(B_seed,np.interp(grid,seed_w,seed_n),rcond=None)[0],lower[:m]+1e-5,upper[:m]-1e-5)
    if settings.initial_thickness_nm is not None and n_thick:
        d_seed=float(settings.initial_thickness_nm)
        if not np.isfinite(d_seed) or d_seed<=0:
            raise ValueError('Invalid Swanepoel thickness seed')
        x0[2*m:]=np.clip(d_seed,lower[2*m:]+1e-5,upper[2*m:]-1e-5)
    def residual(x,include_penalty=True):
        a=x[:m]; b=x[m:2*m]; d=x[2*m:] if n_thick else thicknesses
        out=[]
        for s,B,di in zip(spectra,mats,d):
            n=B@a; k=B@b
            violation_n=np.maximum(1e-6-n,0.)
            violation_k=np.maximum(-k,0.)
            R,T=coherent_film(s.wavelength_nm,np.maximum(n,1e-6)+1j*np.maximum(k,0.),di,s.substrate_n,s.ambient_n,s.angle_deg,s.polarization,roughness_nm=getattr(s,'roughness_nm',0.),substrate_thickness_nm=getattr(s,'substrate_thickness_nm',None),backing_n=getattr(s,'backing_n',1.))
            out.extend(1000.*violation_n)
            out.extend(1000.*violation_k)
            if s.R is not None: out.extend((R-s.R)/s.sigma_R)
            if s.T is not None: out.extend((T-s.T)/s.sigma_T)
        if include_penalty:
            out.extend(np.sqrt(settings.smoothness)*np.diff(a,2))
            out.extend(np.sqrt(settings.smoothness)*np.diff(b,2))
        return np.asarray(out)
    return residual,x0,lower,upper,(lo,hi),thicknesses

def harris_hawks(objective,lower,upper,population=18,iterations=30,seed=42,initial=None):
    '''Bounded HHO: exploration, soft/hard besiege, and rapid dives with Levy steps.'''
    rng=np.random.default_rng(seed); dim=len(lower)
    hawks=rng.uniform(lower,upper,size=(population,dim))
    if initial is not None: hawks[0]=np.clip(np.asarray(initial,dtype=float),lower,upper)
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
        initial,history,evaluations_global=harris_hawks(objective,low,high,settings.population,settings.iterations,settings.seed,initial=x0)
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

def fit_joint_two_stage(spectra, settings):
    """Blind joint Cauchy/log-absorption initialization followed by B-spline refinement."""
    validate_spectra(spectra)
    nominal=np.array([s.thickness_nm for s in spectra],float)
    count=len(spectra) if settings.fit_thickness else 0
    low=np.r_[settings.n_min,-1.0,-1.0,-20.0,-8.0,nominal*(1-settings.thickness_tolerance) if count else []]
    high=np.r_[settings.n_max,1.0,1.0,np.log(settings.k_max),8.0,nominal*(1+settings.thickness_tolerance) if count else []]
    def evaluate(q, penalty=True):
        out=[]
        for j,s in enumerate(spectra):
            x=800.0/np.asarray(s.wavelength_nm,float)
            n=q[0]+q[1]*x*x+q[2]*x**4
            k=np.exp(np.clip(q[3]+q[4]*(x-1),-30,2))
            d=q[5+j] if count else nominal[j]
            # Keep the optimizer objective defined for inadmissible trial dispersions.
            # Penalize the violation rather than silently accepting clipped optics.
            violation=np.maximum(1e-6-n,0.)
            R,T=coherent_film(s.wavelength_nm,np.maximum(n,1e-6)+1j*k,d,s.substrate_n,s.ambient_n,s.angle_deg,s.polarization,roughness_nm=getattr(s,'roughness_nm',0.),substrate_thickness_nm=getattr(s,'substrate_thickness_nm',None),backing_n=getattr(s,'backing_n',1.))
            out.extend(1000.*violation)
            if s.R is not None:out.extend((R-s.R)/s.sigma_R)
            if s.T is not None:out.extend((T-s.T)/s.sigma_T)
        if penalty:out.extend([.001*q[2],.001*q[4]])
        return np.asarray(out)
    seeds=[]
    for n0 in (1.5,1.75,2.0,2.25):
        for b0 in (0.0,.04):
            seeds.append(np.r_[n0,b0,0.,np.log(.001),0.,nominal if count else []])
    trials=[least_squares(evaluate,np.clip(q,low+1e-8,high-1e-8),
             bounds=(low,high),max_nfev=max(120,settings.max_nfev),x_scale='jac') for q in seeds]
    first=min(trials,key=lambda z:np.linalg.norm(evaluate(z.x,False)))
    def param_values(w):
        x=800./np.asarray(w,float)
        return first.x[0]+first.x[1]*x*x+first.x[2]*x**4,np.exp(np.clip(first.x[3]+first.x[4]*(x-1),-30,2))
    residual,x0,lo,hi,domain,_=make_problem(spectra,settings)
    m=settings.n_basis
    grid=np.linspace(*domain,max(4*m,80))
    B=spline_design(grid,*domain,m)
    pn,pk=param_values(grid)
    x0[:m]=np.linalg.lstsq(B,pn,rcond=None)[0]
    x0[m:2*m]=np.linalg.lstsq(B,pk,rcond=None)[0]
    if count:x0[2*m:]=first.x[5:]
    x0=np.clip(x0,lo+1e-8,hi-1e-8)
    second=least_squares(residual,x0,bounds=(lo,hi),max_nfev=max(180,settings.max_nfev),x_scale='jac')
    score1=float(np.sqrt(np.mean(evaluate(first.x,False)**2)))
    score2=float(np.sqrt(np.mean(residual(second.x,False)**2)))
    use_spline=score2<=score1*(1+1e-4)
    chosen=second if use_spline else first
    def optical_constants(w):
        if not use_spline:return param_values(w)
        b=spline_design(np.asarray(w,float),*domain,m)
        return b@second.x[:m],b@second.x[m:2*m]
    thickness=(second.x[2*m:] if use_spline else first.x[5:]) if count else nominal
    return {'result':chosen,'optical_constants':optical_constants,'thicknesses':thickness,
            'weighted_rmse':score2 if use_spline else score1,'history':[],
            'global_evaluations':sum(z.nfev for z in trials),'domain':domain,
            'stages':{'parametric_weighted_rmse':score1,'spline_weighted_rmse':score2,
                      'selected':'spline' if use_spline else 'parametric',
                      'parametric_thicknesses_nm':(first.x[5:] if count else nominal).tolist(),
                      'spline_thicknesses_nm':(second.x[2*m:] if count else nominal).tolist(),
                      'parametric_coefficients':dict(a=float(first.x[0]),b=float(first.x[1]),c=float(first.x[2]),
                                                     q0=float(first.x[3]),q1=float(first.x[4]))}}

def simulate_four(wavelength_nm,n,k):
    specs=[]
    for name,d,sub in [('A thick air',1000.,1.),('B thin air',50.,1.),('C thick glass',1000.,1.5),('D thin glass',50.,1.5)]:
        R,T=coherent_film(wavelength_nm,n+1j*k,d,sub)
        specs.append(Spectrum(np.array(wavelength_nm),R,T,d,substrate_n=sub,name=name))
    return specs

def compare_prediction(spectrum,n,k,thickness=None):
    R,T=coherent_film(spectrum.wavelength_nm,np.asarray(n)+1j*np.asarray(k),spectrum.thickness_nm if thickness is None else thickness,spectrum.substrate_n,spectrum.ambient_n,spectrum.angle_deg,spectrum.polarization,roughness_nm=getattr(spectrum,'roughness_nm',0.),substrate_thickness_nm=getattr(spectrum,'substrate_thickness_nm',None),backing_n=getattr(spectrum,'backing_n',1.))
    diffs=np.concatenate([R-spectrum.R if spectrum.R is not None else np.array([]),T-spectrum.T if spectrum.T is not None else np.array([])])
    return {'R':R,'T':T,'rmse':float(np.sqrt(np.mean(diffs**2)))}

def phase2_hybrid_seed(wavelength_nm, global_wavelength_nm, global_n, global_k, phase1_parameters, nodes=25, transition_nm=690.0):
    """Reproduce the original phase2_fair_fit.py knot initialization.

    Below and at 690 nm interpolate the independent broadband global estimate.
    Above 690 nm use Phase I with a fixed 800 nm reference wavelength.
    The broadband input must come from an independently calculated fit, not
    from a reference n/k table or an unrelated material.
    """
    w=np.asarray(wavelength_nm,dtype=float)
    gw=np.asarray(global_wavelength_nm,dtype=float)
    gn=np.asarray(global_n,dtype=float)
    gk=np.asarray(global_k,dtype=float)
    if w.ndim!=1 or gw.ndim!=1 or len(w)<nodes or len(gw)<2 or np.any(np.diff(w)<=0) or np.any(np.diff(gw)<=0):
        raise ValueError('Hybrid seed requires increasing wavelength arrays')
    if not (np.all(np.isfinite(gn)) and np.all(np.isfinite(gk)) and len(gw)==len(gn)==len(gk)) or np.any(gk<=0):
        raise ValueError('Broadband seed requires finite n and positive k')
    kn=np.linspace(w[0],w[-1],nodes)
    if gw[0]>kn[0] or gw[-1]<kn[-1]:
        raise ValueError('Broadband seed must cover the complete Phase II wavelength range')
    n0=np.interp(kn,gw,gn)
    k0=np.interp(kn,gw,gk)
    for j,wl in enumerate(kn):
        if wl>transition_nm:
            x=800.0/wl
            p=phase1_parameters
            n0[j]=p['a']+p['b']*x*x+p['c']*x**4
            k0[j]=np.exp(np.clip(p['q0']+p['q1']*(x-1),-35,3))
    return np.r_[n0,np.log(np.maximum(k0,1e-8))]


def estimate_broadband_seed(spectrum, n_seed, k_seed, nodes=25, max_nfev=55):
    """Independent broadband nodal fit from the uploaded R/T measurements.

    No material-specific reference curves are used. A smooth, bounded fit
    supplies the low-wavelength portion of the historical hybrid initializer.
    """
    from scipy.interpolate import CubicSpline
    w=np.asarray(spectrum.wavelength_nm,dtype=float)
    if spectrum.R is None or spectrum.T is None:
        raise ValueError('Broadband initialization requires measured R and T')
    kn=np.linspace(w[0],w[-1],nodes)
    n0=np.interp(kn,w,np.asarray(n_seed,dtype=float))
    k0=np.interp(kn,w,np.asarray(k_seed,dtype=float))
    q0=np.r_[np.clip(n0,1.01,3.99),np.log(np.clip(k0,1e-7,2.0))]
    lo=np.r_[np.full(nodes,1.),np.full(nodes,-20.)]
    hi=np.r_[np.full(nodes,4.),np.full(nodes,1.)]
    def objective(q):
        n=np.maximum(CubicSpline(kn,q[:nodes],bc_type='natural')(w),1e-6)
        k=np.exp(np.clip(CubicSpline(kn,q[nodes:],bc_type='natural')(w),-24,2))
        R,T=coherent_film(w,n+1j*k,spectrum.thickness_nm,spectrum.substrate_n,spectrum.ambient_n,spectrum.angle_deg,spectrum.polarization,roughness_nm=getattr(spectrum,'roughness_nm',0.),substrate_thickness_nm=getattr(spectrum,'substrate_thickness_nm',None),backing_n=getattr(spectrum,'backing_n',1.))
        return np.r_[T-spectrum.T,R-spectrum.R,.012*np.diff(q[:nodes],2),.003*np.diff(q[nodes:],2)]
    sol=least_squares(objective,np.clip(q0,lo+1e-8,hi-1e-8),bounds=(lo,hi),max_nfev=max_nfev,ftol=1e-6,xtol=1e-6,gtol=1e-6)
    n=np.maximum(CubicSpline(kn,sol.x[:nodes],bc_type='natural')(w),1e-6)
    k=np.exp(np.clip(CubicSpline(kn,sol.x[nodes:],bc_type='natural')(w),-24,2))
    return {'n':n,'k':k,'nfev':int(sol.nfev),'success':bool(sol.success),'message':str(sol.message)}


def fit_historical_logk(spectrum, max_nfev=90, nodes=25, n_seed=None, k_seed=None, warm_start=None):
    """Phase-II natural cubic nodal splines: n and log(k), fixed known thickness.
    Recovered phase2_fair_fit.py uses 25 nodes, residual penalties .009/.0015.
    This implementation does NOT reproduce its archived warm-start unless supplied.
    """
    from scipy.interpolate import CubicSpline
    validate_spectra([spectrum])
    w=np.asarray(spectrum.wavelength_nm,dtype=float)
    if len(w)<nodes+2:
        raise ValueError('Historical spline needs at least 27 wavelength samples')
    if spectrum.R is None or spectrum.T is None:
        raise ValueError('Historical spline requires both measured R and T')
    kn=np.linspace(w[0],w[-1],nodes)
    def values(q, wl):
        nn=np.maximum(CubicSpline(kn,q[:nodes],bc_type='natural')(wl),1e-6)
        kk=np.exp(np.clip(CubicSpline(kn,q[nodes:],bc_type='natural')(wl),-24,2))
        return nn,kk
    initial_n=np.full(nodes,1.75) if n_seed is None else np.interp(kn,w,np.asarray(n_seed))
    initial_k=np.full(nodes,.05) if k_seed is None else np.interp(kn,w,np.asarray(k_seed))
    x0=np.r_[np.clip(initial_n,1.001,3.999),np.log(np.clip(initial_k,np.exp(-19.9),2.5))]
    if warm_start is not None:
        q=np.asarray(warm_start,dtype=float)
        if q.shape != (2*nodes,) or not np.all(np.isfinite(q)):
            raise ValueError('Warm start must contain 2*nodes finite spline coefficients')
        x0=q.copy()
    lower=np.r_[np.full(nodes,1.),np.full(nodes,-20.)]
    upper=np.r_[np.full(nodes,4.),np.full(nodes,1.)]
    train=np.arange(len(w))%5!=0
    holdout=~train
    def predict(q):
        n,k=values(q,w)
        return coherent_film(w,n+1j*k,spectrum.thickness_nm,spectrum.substrate_n,spectrum.ambient_n,spectrum.angle_deg,spectrum.polarization,roughness_nm=getattr(spectrum,'roughness_nm',0.),substrate_thickness_nm=getattr(spectrum,'substrate_thickness_nm',None),backing_n=getattr(spectrum,'backing_n',1.))
    def residual(q):
        R,T=predict(q)
        return np.r_[(T-spectrum.T)[train],(R-spectrum.R)[train],.009*np.diff(q[:nodes],2),.0015*np.diff(q[nodes:],2)]
    if max_nfev == 0:
        from types import SimpleNamespace
        result=SimpleNamespace(x=x0,success=True,message='Archived coefficients evaluated without optimization',nfev=0)
    else:
        # Distinct physically admissible starts mitigate local-minimum dependence.
        # All candidates use the identical train/holdout split and regularization.
        seeds=[x0]
        if warm_start is not None:
            seeds.append(np.r_[np.full(nodes,1.75),np.full(nodes,np.log(.03))])
        else:
            seeds.extend([np.r_[np.full(nodes,1.7),np.full(nodes,np.log(.03))],
                          np.r_[np.full(nodes,1.9),np.full(nodes,np.log(.12))]])
        trials=[]
        for seed in seeds:
            candidate=least_squares(residual,np.clip(seed,lower+1e-9,upper-1e-9),
                bounds=(lower,upper),max_nfev=max_nfev,ftol=1e-7,xtol=1e-7,gtol=1e-7)
            trials.append(candidate)
        result=min(trials,key=lambda trial: float(np.sum(residual(trial.x)**2)))

    n,k=values(result.x,w)
    R,T=predict(result.x)
    def rmse(mask):
        return float(np.sqrt(np.mean(np.r_[(T-spectrum.T)[mask],(R-spectrum.R)[mask]]**2)))
    return {'result':result,'n':n,'k':k,'R':R,'T':T,
            'train_rmse':rmse(train),'holdout_rmse':rmse(holdout),
            'train_mask':train,'nodes':nodes,'initialization':('saved coefficients' if warm_start is not None else 'provided n/k seed' if n_seed is not None else 'flat seed'),'success':bool(result.success),'message':str(result.message),'nfev':int(result.nfev),'candidate_count':(len(trials) if max_nfev!=0 else 1),'regularized_objective':float(np.sum(residual(result.x)**2))}


def fit_absorption_aware_blind(spectrum, settings=None, transparent_min_nm=None, starts=5):
    """Data-only Cauchy/attenuation initialization followed by bounded spline R/T inversion.

    A small deterministic set of index-dispersion and absorption seeds is ranked
    using the FULL R/T objective, not just the weak-absorption window.
    No known optical constants enter the initializer.
    """
    from scipy.ndimage import gaussian_filter1d
    if spectrum.R is None or spectrum.T is None:
        raise ValueError('Blind absorption-aware initialization needs both R and T')
    settings = settings or FitSettings(n_basis=18, smoothness=.01, max_nfev=140)
    if settings.fit_thickness:
        raise ValueError('Blind initializer requires fixed thickness')
    residual, _, low, high, domain, _ = make_problem([spectrum], settings)
    w=np.asarray(spectrum.wavelength_nm,float)
    R=np.asarray(spectrum.R,float)
    T=np.asarray(spectrum.T,float)
    m=settings.n_basis
    B=spline_design(w,*domain,m)
    threshold=float(transparent_min_nm) if transparent_min_nm is not None else float(np.quantile(w,.35))
    mask=(w>=threshold)&(T>.15)
    if np.count_nonzero(mask)<12:
        raise ValueError('Insufficient weakly absorbing wavelengths for fringe-based initialization')
    substrate=np.asarray(spectrum.substrate_n)
    ns_mask=substrate[mask] if substrate.ndim else spectrum.substrate_n
    # Beer-Lambert estimate neglects interfaces; it is only a starting point.
    k_beer=gaussian_filter1d(np.clip(-w*np.log(np.clip(T,1e-5,.999))/
        (4*np.pi*spectrum.thickness_nm),1e-6,1.5),2)
    seeds=[]
    # Fit both Cauchy dispersion and a constant nonnegative k in the selected
    # window. Allowing weak absorption avoids forcing a transparent model onto
    # an only approximately transparent measured region.
    for n0 in np.arange(1.35,2.66,.13):
        def fringe(q):
            n=q[0]+q[1]*(800./w[mask])**2
            rr,tt=coherent_film(w[mask],n+1j*q[2],spectrum.thickness_nm,
                ns_mask,spectrum.ambient_n,spectrum.angle_deg,spectrum.polarization)
            return np.r_[rr-R[mask],tt-T[mask]]
        trial=least_squares(fringe,[n0,0.,.01],bounds=([1.,-.9,0.],[3.8,.9,.4]),
            max_nfev=75)
        seeds.append(trial.x)
    # Retain diverse seeds across n values to avoid the window's local fringe aliases.
    candidates=[]
    for q in seeds:
        n_est=q[0]+q[1]*(800./w)**2
        for k_factor in (.5,1.,1.5):
            k_est=np.clip(k_factor*k_beer,1e-6,1.8)
            a=np.linalg.lstsq(B,n_est,rcond=None)[0]
            b=np.linalg.lstsq(B,k_est,rcond=None)[0]
            x=np.clip(np.r_[a,b],low+1e-7,high-1e-7)
            candidates.append((float(np.dot(residual(x),residual(x))),x))
    candidates.sort(key=lambda item:item[0])
    chosen=[]
    for _,x in candidates:
        if not any(np.linalg.norm(x[:m]-other[:m])<.05 for other in chosen):
            chosen.append(x)
        if len(chosen)>=max(1,int(starts)):
            break
    fits=[]
    for x in chosen:
        fit=least_squares(residual,x,bounds=(low,high),
            max_nfev=settings.max_nfev,x_scale='jac')
        fits.append(fit)
    result=min(fits,key=lambda fit:float(np.dot(fit.fun,fit.fun)))
    n=B@result.x[:m];k=B@result.x[m:2*m]
    rp,tp=coherent_film(w,n+1j*k,spectrum.thickness_nm,spectrum.substrate_n,
        spectrum.ambient_n,spectrum.angle_deg,spectrum.polarization)
    return dict(n=n,k=k,R_fit=rp,T_fit=tp,thickness_nm=spectrum.thickness_nm,
        R_rmse=float(np.sqrt(np.mean((rp-R)**2))),
        T_rmse=float(np.sqrt(np.mean((tp-T)**2))),
        nfev=int(sum(fit.nfev for fit in fits)),success=bool(result.success),
        candidates=len(fits),transparent_threshold_nm=threshold)
