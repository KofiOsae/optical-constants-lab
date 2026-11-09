"""Normal-incidence coherent TMM for an effective homogeneous film.

This calculates specular R and T only. It does not model diffuse scattering.
Time convention exp(-i omega t), complex n+ik, k>=0.
"""
import numpy as np

def film_rt(wavelength_nm,film_nk,thickness_nm,substrate_n=1.5,incident_n=1.0):
    wl=np.asarray(wavelength_nm,float)
    nf=np.asarray(film_nk,complex)
    ns=np.asarray(substrate_n,complex)
    n0=complex(incident_n)
    wl,nf,ns=np.broadcast_arrays(wl,nf,ns)
    d=float(thickness_nm)
    if (not np.isfinite(d) or d<0 or np.any(~np.isfinite(wl)) or np.any(wl<=0)
        or np.any(~np.isfinite(nf)) or np.any(nf.imag<0)
        or np.any(~np.isfinite(ns)) or np.any(ns.imag<0)
        or not np.isfinite(n0) or n0.imag!=0 or n0.real<=0 or np.any(ns.real<=0)):
        raise ValueError("Invalid passive film, substrate or thickness")
    r01=(n0-nf)/(n0+nf)
    r12=(nf-ns)/(nf+ns)
    t01=2*n0/(n0+nf)
    t12=2*nf/(nf+ns)
    phase=np.exp(2j*np.pi*nf*d/wl)
    denom=1+r01*r12*phase**2
    r=(r01+r12*phase**2)/denom
    t=t01*t12*phase/denom
    R=np.abs(r)**2
    T=np.real(ns)/n0.real*np.abs(t)**2
    A=1-R-T
    if np.any(~np.isfinite(R+T)) or np.any(A<-1e-7):
        raise ValueError("Nonphysical TMM energy balance")
    return {"R":R,"T":T,"A":np.maximum(0,A)}
