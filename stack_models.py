"""Coherent stratified-media Fresnel solver and ellipsometric observables.
Convention exp(-i omega t), passive refractive index n+ik. Semi-infinite
incident and exit media must be lossless for absolute power T.
"""
import numpy as np

def stack_amplitudes(wavelength_nm, layers, ambient_n=1., substrate_n=1.5, angle_deg=0., polarization="s"):
    wl=np.asarray(wavelength_nm,float)
    if wl.ndim!=1 or np.any(wl<=0): raise ValueError("Wavelengths must be positive 1D nm")
    ns=[complex(ambient_n)]+[complex(l["n"])+1j*float(l.get("k",0)) for l in layers]+[complex(substrate_n)]
    if abs(ns[0].imag)>1e-12 or abs(ns[-1].imag)>1e-12: raise ValueError("External media must be lossless")
    if any(z.real<=0 or z.imag<0 for z in ns): raise ValueError("Require n>0, k>=0")
    if polarization not in ("s","p"): raise ValueError("Choose s or p")
    u=ns[0]*np.sin(np.deg2rad(float(angle_deg)))
    q=[np.sqrt(z*z-u*u+0j) for z in ns]
    q=[(-v if v.imag<0 else v) for v in q]
    y=q if polarization=="s" else [z*z/v for z,v in zip(ns,q)]
    # Recursive effective reflection/transmission from substrate upwards.
    r=np.zeros(wl.shape,complex);t=np.ones(wl.shape,complex)
    for j in range(len(ns)-2,-1,-1):
        rij=(y[j]-y[j+1])/(y[j]+y[j+1])
        tij=2*y[j]/(y[j]+y[j+1])
        if j+1<=len(layers):
            thickness=float(layers[j]["thickness_nm"])
            if thickness<0: raise ValueError("Layer thickness must be nonnegative")
            phase=np.exp(2j*np.pi*q[j+1]*thickness/wl)
        else: phase=np.ones_like(wl,complex)
        denominator=1+rij*r*phase*phase
        t=tij*phase*t/denominator
        r=(rij+r*phase*phase)/denominator
    return r,t,float(np.real(y[-1]/y[0]))

def stack_spectrum(wavelength_nm,layers,ambient_n=1.,substrate_n=1.5,angle_deg=0.,polarization="s"):
    r,t,f=stack_amplitudes(wavelength_nm,layers,ambient_n,substrate_n,angle_deg,polarization)
    return {"wavelength_nm":list(map(float,wavelength_nm)),"R":np.abs(r)**2,"T":f*np.abs(t)**2,"A":1-np.abs(r)**2-f*np.abs(t)**2}

def ellipsometry(wavelength_nm,layers,ambient_n=1.,substrate_n=1.5,angle_deg=70.):
    rp,_,_=stack_amplitudes(wavelength_nm,layers,ambient_n,substrate_n,angle_deg,"p")
    rs,_,_=stack_amplitudes(wavelength_nm,layers,ambient_n,substrate_n,angle_deg,"s")
    rho=rp/np.where(np.abs(rs)>1e-12,rs,np.nan+0j)
    return {"wavelength_nm":list(map(float,wavelength_nm)),"psi_deg":np.rad2deg(np.arctan(np.abs(rho))),"delta_deg":np.rad2deg(np.angle(rho)),"rho_real":rho.real,"rho_imag":rho.imag}
