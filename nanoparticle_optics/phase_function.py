"""Polarized far-field angular scattering of homogeneous Mie spheres."""
import numpy as np
from nanoparticle_optics.angular import _coefficients

def sphere_phase_function(wavelength_nm,radius_nm,particle_index,host_index=1.,
                          angles_deg=None):
    """Return normalized phase function p(theta), integral p dOmega = 1.

    S1 is perpendicular and S2 parallel to the scattering plane.
    Phase functions are differential scattering cross sections normalized
    to total scattering; an unpolarized incident plane wave is assumed.
    """
    theta=np.linspace(0,180,181) if angles_deg is None else np.asarray(angles_deg,float)
    if theta.ndim!=1 or len(theta)<3 or not np.all(np.isfinite(theta)) or np.any(np.diff(theta)<=0) or theta[0]<0 or theta[-1]>180:
        raise ValueError("Angles must be a strictly increasing grid within 0–180 degrees")
    x,ell,a,b=_coefficients(wavelength_nm,radius_nm,particle_index,host_index)
    mu=np.cos(np.deg2rad(theta))
    pi0=np.zeros_like(mu);pi1=np.ones_like(mu)
    S1=np.zeros_like(mu,dtype=complex);S2=np.zeros_like(mu,dtype=complex)
    for j in range(1,len(ell)+1):
        tau=j*mu*pi1-(j+1)*pi0
        factor=(2*j+1)/(j*(j+1))
        S1+=factor*(a[j-1]*pi1+b[j-1]*tau)
        S2+=factor*(a[j-1]*tau+b[j-1]*pi1)
        pi2=((2*j+1)*mu*pi1-(j+1)*pi0)/j
        pi0,pi1=pi1,pi2
    # dC/dOmega = (|S1|^2+|S2|^2)/(2 k_host^2)
    k=2*np.pi*host_index/float(wavelength_nm)
    perp=np.abs(S1)**2/(k*k)
    para=np.abs(S2)**2/(k*k)
    total=.5*(perp+para)
    csca=2*np.pi*np.trapz(total*np.sin(np.deg2rad(theta)),np.deg2rad(theta))
    if csca<=0 or not np.isfinite(csca):raise ValueError("Angular scattering normalization undefined")
    phase=total/csca
    return {"angle_deg":theta,"phase_sr_inv":phase,
            "perpendicular_nm2_sr":perp,"parallel_nm2_sr":para,
            "Csca_angular_nm2":float(csca),
            "g_angular":float(2*np.pi*np.trapz(phase*mu*np.sin(np.deg2rad(theta)),np.deg2rad(theta)))}
