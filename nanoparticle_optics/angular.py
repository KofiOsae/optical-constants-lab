"""Angular Mie scattering diagnostics for homogeneous spheres in a lossless host.

Uses the exp(-i omega t) convention of nanoparticle_optics.mie.
The anisotropy is obtained from the exact Mie-coefficient sum, rather
than a Rayleigh approximation or a numerical derivative of total Qsca.
"""
import numpy as np
from scipy.special import spherical_jn, spherical_yn
from nanoparticle_optics.mie import mie_sphere


def _coefficients(wavelength_nm, radius_nm, particle_index, host_index):
    wl=float(wavelength_nm);radius=float(radius_nm)
    n=complex(particle_index);nh=float(host_index)
    if not all(np.isfinite(v) for v in (wl,radius,n.real,n.imag,nh)) or wl<=0 or radius<=0 or nh<=0 or n.imag<0:
        raise ValueError("Invalid sphere, wavelength or optical constants")
    x=2*np.pi*nh*radius/wl;m=n/nh
    if x>1000 or abs(m*x)>1000:raise ValueError("Size parameter outside validated numerical range")
    count=max(3,int(np.ceil(x+4*x**(1/3)+2)))
    ell=np.arange(1,count+1)
    z=m*x
    jx=spherical_jn(ell,x);jxp=spherical_jn(ell,x,derivative=True)
    yx=spherical_yn(ell,x);yxp=spherical_yn(ell,x,derivative=True)
    jz=spherical_jn(ell,z);jzp=spherical_jn(ell,z,derivative=True)
    px=x*jx;pxp=jx+x*jxp
    pz=z*jz;pzp=jz+z*jzp
    xi=x*(jx+1j*yx);xip=(jx+1j*yx)+x*(jxp+1j*yxp)
    with np.errstate(divide="ignore",invalid="ignore",over="ignore"):
        a=(m*pz*pxp-px*pzp)/(m*pz*xip-xi*pzp)
        b=(pz*pxp-m*px*pzp)/(pz*xip-m*xi*pzp)
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("Mie coefficient numerical instability")
    return x,ell,a,b


def mie_anisotropy(wavelength_nm,radius_nm,particle_index,host_index=1.0):
    """Return Qsca, g and Csca for a single sphere at each wavelength.

    g=0 is used when Csca=0 (direction undefined); the flag identifies it.
    """
    wl=np.asarray(wavelength_nm,float)
    index=np.asarray(particle_index,complex)
    wl,index=np.broadcast_arrays(wl,index)
    out=np.empty(wl.shape);qs=np.empty(wl.shape)
    undefined=np.zeros(wl.shape,dtype=bool)
    for idx in np.ndindex(wl.shape):
        x,ell,a,b=_coefficients(wl[idx],radius_nm,index[idx],host_index)
        q=2/x**2*np.sum((2*ell+1)*(np.abs(a)**2+np.abs(b)**2))
        qs[idx]=q
        if q<=1e-24:
            out[idx]=0.;undefined[idx]=True;continue
        coupling=np.sum(ell[:-1]*(ell[:-1]+2)/(ell[:-1]+1)*
                       np.real(a[:-1]*np.conj(a[1:])+b[:-1]*np.conj(b[1:])))
        cross=np.sum((2*ell+1)/(ell*(ell+1))*np.real(a*np.conj(b)))
        out[idx]=4*(coupling+cross)/(x*x*q)
        if not np.isfinite(out[idx]) or out[idx]<-1-1e-8 or out[idx]>1+1e-8:
            raise ValueError("Nonphysical Mie anisotropy; numerical convergence failed")
        out[idx]=np.clip(out[idx],-1.,1.)
    return {"wavelength_nm":wl,"g":out,"Qsca":qs,
            "Csca_nm2":qs*np.pi*float(radius_nm)**2,"g_undefined":undefined}


def ensemble_anisotropy(wavelength_nm,radii_nm,number_weights,particle_index,host_index=1.0,
                        number_density_nm3=None):
    """Number-average Csca and scattering-weighted g for independent spheres.

    mu_s = N*<Csca>, mu_s_prime = N*<Csca*(1-g)> in inverse nm.
    """
    wl=np.asarray(wavelength_nm,float)
    r=np.asarray(radii_nm,float);w=np.asarray(number_weights,float)
    if wl.ndim!=1 or len(wl)==0 or np.any(~np.isfinite(wl)) or np.any(wl<=0):
        raise ValueError("Invalid wavelength grid")
    if r.ndim!=1 or r.size==0 or w.shape!=r.shape or np.any(~np.isfinite(r)) or np.any(r<=0) or np.any(~np.isfinite(w)) or np.any(w<0) or w.sum()<=0:
        raise ValueError("Invalid radii or number weights")
    w=w/w.sum()
    sca=np.zeros_like(wl);moment=np.zeros_like(wl)
    for radius,weight in zip(r,w):
        if weight==0:continue
        v=mie_anisotropy(wl,float(radius),particle_index,host_index)
        c=v["Csca_nm2"]
        sca+=weight*c
        moment+=weight*c*v["g"]
    g=np.divide(moment,sca,out=np.zeros_like(sca),where=sca>0)
    result={"wavelength_nm":wl,"g":g,"Csca_nm2":sca,
            "Csca_transport_nm2":np.maximum(0,sca-moment),
            "g_undefined":sca<=0}
    if number_density_nm3 is not None:
        density=float(number_density_nm3)
        if not np.isfinite(density) or density<0:raise ValueError("Number density must be nonnegative")
        result["mu_s_nm_inv"]=density*sca
        result["mu_s_prime_nm_inv"]=density*np.maximum(0,sca-moment)
    return result
