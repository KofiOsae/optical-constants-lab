"""Coherent planar multilayer Fresnel recursion (exp(-i omega t), N=n+ik).
Incident and exit half-spaces must be lossless for power transmittance.
All layers are coherent, isotropic, homogeneous and planar.
"""
import numpy as np

def multilayer_rt(wavelength_nm, layers, substrate_n=1.46, ambient_n=1.0,
                  angle_deg=0.0, polarization='s'):
    wl=np.atleast_1d(np.asarray(wavelength_nm,float))
    if np.any(~np.isfinite(wl)) or np.any(wl<=0):raise ValueError('Wavelengths must be positive and finite')
    if polarization not in ('s','p'):raise ValueError('polarization must be s or p')
    def expand(n):
        a=np.asarray(n,complex)
        return np.broadcast_to(a,wl.shape)
    n0=expand(ambient_n);ns=expand(substrate_n)
    if np.any(abs(n0.imag)>1e-12) or np.any(abs(ns.imag)>1e-12):
        raise ValueError('Power transmittance currently requires lossless ambient and exit substrate')
    media=[n0]
    ds=[]
    for layer in layers:
        n=expand(layer['n']);d=float(layer['thickness_nm'])
        if d<0 or not np.isfinite(d):raise ValueError('Invalid layer thickness')
        if np.any(n.real<=0) or np.any(n.imag<0):raise ValueError('Require passive layer n>0,k>=0')
        media.append(n);ds.append(d)
    media.append(ns)
    u=n0*np.sin(np.deg2rad(angle_deg))
    qs=[]
    for n in media:
        q=np.sqrt(n*n-u*u+0j)
        q=np.where((q.imag<0)|((abs(q.imag)<1e-13)&(q.real<0)),-q,q)
        qs.append(q)
    ys=[q if polarization=='s' else n*n/q for n,q in zip(media,qs)]
    r=[];t=[]
    for i in range(len(media)-1):
        r.append((ys[i]-ys[i+1])/(ys[i]+ys[i+1]))
        t.append(2*ys[i]/(ys[i]+ys[i+1]))
    re=r[-1];te=t[-1]
    for i in range(len(ds)-1,-1,-1):
        phase=np.exp(1j*2*np.pi*qs[i+1]*ds[i]/wl)
        den=1+r[i]*re*phase*phase
        te=t[i]*te*phase/den
        re=(r[i]+re*phase*phase)/den
    R=np.abs(re)**2
    T=np.real(ys[-1]/ys[0])*np.abs(te)**2
    return R,T
