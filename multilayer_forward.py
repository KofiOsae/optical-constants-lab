"""Coherent planar-stack R/T kernel for future inverse characterization.

Layer order is incident ambient -> finite layers -> semi-infinite substrate.
Conventions: exp(-i omega t), N=n+i*k, positive forward attenuation.
This is a minimal inversion forward kernel, not a coating-design interface.
"""
import numpy as np

def coherent_stack(wavelength_nm, layers, substrate_n=1.5, ambient_n=1.,
                   angle_deg=0., polarization='s'):
    """Return power R,T for isotropic coherent finite layers.

    layers: sequence of (complex refractive index, thickness_nm).
    Each index may be scalar or wavelength-length array. Ambient and
    substrate must be lossless; finite layers must be passive.
    """
    w=np.asarray(wavelength_nm,dtype=float)
    if w.ndim!=1 or w.size==0 or not np.all(np.isfinite(w)) or np.any(w<=0):
        raise ValueError('Wavelengths must be a nonempty positive finite 1-D array')
    def optical(value):
        a=np.asarray(value,dtype=complex)
        if a.ndim>1 or (a.ndim==1 and a.shape!=w.shape):
            raise ValueError('Optical index must be scalar or wavelength-aligned')
        return np.broadcast_to(a,w.shape)
    n0=optical(ambient_n);ns=optical(substrate_n)
    if np.any(np.abs(n0.imag)>1e-12) or np.any(np.abs(ns.imag)>1e-12):
        raise ValueError('Ambient and substrate must be lossless')
    if np.any(n0.real<=0) or np.any(ns.real<=0):
        raise ValueError('Ambient and substrate indices must be positive')
    finite=[]
    for index,d in layers:
        n=optical(index)
        if not np.isfinite(d) or d<0 or np.any(~np.isfinite(n)) or np.any(n.real<=0) or np.any(n.imag<0):
            raise ValueError('Finite layer needs passive index and nonnegative thickness')
        finite.append((n,float(d)))
    n=[n0]+[v for v,_ in finite]+[ns]
    u=n0*np.sin(np.deg2rad(angle_deg))
    q=[np.sqrt(v*v-u*u+0j) for v in n]
    q=[np.where((x.imag<0)|((np.abs(x.imag)<1e-14)&(x.real<0)),-x,x) for x in q]
    pol=polarization.lower()
    if pol=='s':y=q
    elif pol=='p':y=[v*v/x for v,x in zip(n,q)]
    else:raise ValueError('polarization must be s or p')
    reflection=[];transmission=[]
    for j in range(len(n)-1):
        reflection.append((y[j]-y[j+1])/(y[j]+y[j+1]))
        transmission.append(2*y[j]/(y[j]+y[j+1]))
    r=reflection[-1];t=transmission[-1]
    for j in range(len(n)-3,-1,-1):
        phase=np.exp(2j*np.pi*q[j+1]*finite[j][1]/w)
        den=1+reflection[j]*r*phase*phase
        t=transmission[j]*t*phase/den
        r=(reflection[j]+r*phase*phase)/den
    R=np.abs(r)**2
    T=np.real(y[-1]/y[0])*np.abs(t)**2
    return R,T
