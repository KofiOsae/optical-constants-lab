"""Three-component quasistatic EMA for connected host, spherical particles and pores.

Sequential Maxwell-Garnett: particles mixed into host first, then pores into
that effective matrix. Order dependent; not a symmetric three-phase solution.
For Bruggeman, solve three-phase implicit relation using continuation.
"""
import numpy as np
from scipy.optimize import root
from nanoparticle_optics.effective_medium import maxwell_garnett

def three_phase_epsilon(host,particle,pore,f_particle,f_pore,method="maxwell_garnett"):
    fp=float(f_particle);fv=float(f_pore)
    if not np.isfinite(fp+fv) or fp<0 or fv<0 or fp+fv>1:
        raise ValueError("Particle and pore fractions must be nonnegative and sum to at most 1")
    h,p,v=np.broadcast_arrays(np.asarray(host,complex),np.asarray(particle,complex),np.asarray(pore,complex))
    if np.any(~np.isfinite(h+p+v)) or np.any(np.imag(h)<-1e-12) or np.any(np.imag(p)<-1e-12) or np.any(np.imag(v)<-1e-12):
        raise ValueError("Finite passive constituent permittivities required")
    if method=="maxwell_garnett":
        if fp+fv==1 and fp>0 and fv>0:
            raise ValueError("Sequential Maxwell-Garnett needs a nonzero connected host fraction")
        matrix=maxwell_garnett(h,p,fp/(1-fv)) if fv<1 else h
        result=maxwell_garnett(matrix,v,fv)
    elif method=="bruggeman":
        result=np.empty(h.shape,complex)
        flat=[a.ravel() for a in (h,p,v)]
        previous=None
        for i,(eh,ep,ev) in enumerate(zip(*flat)):
            if fp==0 and fv==0:sol=eh
            else:
                def fun(x):
                    e=complex(x[0],x[1])
                    residual=(1-fp-fv)*(eh-e)/(eh+2*e)+fp*(ep-e)/(ep+2*e)+fv*(ev-e)/(ev+2*e)
                    return [residual.real,residual.imag]
                guess=previous if previous is not None else (1-fp-fv)*eh+fp*ep+fv*ev
                fit=root(fun,[guess.real,guess.imag],tol=1e-11)
                if not fit.success or np.linalg.norm(fun(fit.x))>1e-7:
                    raise ValueError("Three-phase Bruggeman root did not converge")
                sol=complex(*fit.x)
                if sol.imag < -1e-8:raise ValueError("Nonpassive Bruggeman root")
            result.ravel()[i]=sol
            previous=sol
    else:raise ValueError("Unknown three-phase EMA model")
    if np.any(~np.isfinite(result)) or np.any(result.imag < -1e-8):
        raise ValueError("Nonphysical three-phase effective permittivity")
    return result.item() if result.ndim==0 else result
