"""Homogeneous, nonmagnetic sphere Mie scattering in a lossless host.

Convention exp(-i omega t), N=n+ik with k>=0. Wavelength and radius
are in the same units. The host refractive index must be real and >0.
Outputs are geometric-area-normalized efficiencies; no dependent
scattering or substrate is included. Numerical convergence must be checked
for strongly absorbing or very large spheres.
"""
import numpy as np
from scipy.special import spherical_jn, spherical_yn


def mie_sphere(wavelength_nm, radius_nm, n_particle, n_host=1.0):
    wl, radius, npart, nhost = np.broadcast_arrays(
        np.asarray(wavelength_nm, float), np.asarray(radius_nm, float),
        np.asarray(n_particle, complex), np.asarray(n_host, float))
    if (np.any(~np.isfinite(wl)) or np.any(wl <= 0)
        or np.any(~np.isfinite(radius)) or np.any(radius <= 0)
        or np.any(~np.isfinite(npart)) or np.any(npart.imag < 0)
        or np.any(~np.isfinite(nhost)) or np.any(nhost <= 0)):
        raise ValueError("Invalid wavelength, radius, passive particle index or host index")
    qext = np.empty(wl.shape); qsca = np.empty(wl.shape)
    for idx in np.ndindex(wl.shape):
        x = 2*np.pi*nhost[idx]*radius[idx]/wl[idx]
        m = npart[idx]/nhost[idx]
        if x > 1000 or abs(m*x) > 1000:
            raise ValueError("Size parameter outside validated numerical range")
        count = max(3, int(np.ceil(x + 4*x**(1/3) + 2)))
        ell = np.arange(1, count+1)
        z = m*x
        jx = spherical_jn(ell, x)
        jx_p = spherical_jn(ell, x, derivative=True)
        yx = spherical_yn(ell, x)
        yx_p = spherical_yn(ell, x, derivative=True)
        jz = spherical_jn(ell, z)
        jz_p = spherical_jn(ell, z, derivative=True)
        psi_x = x*jx; psi_x_p = jx+x*jx_p
        psi_z = z*jz; psi_z_p = jz+z*jz_p
        xi_x = x*(jx+1j*yx)
        xi_x_p = (jx+1j*yx)+x*(jx_p+1j*yx_p)
        with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
            a = (m*psi_z*psi_x_p-psi_x*psi_z_p)/(m*psi_z*xi_x_p-xi_x*psi_z_p)
            b = (psi_z*psi_x_p-m*psi_x*psi_z_p)/(psi_z*xi_x_p-m*xi_x*psi_z_p)
        weight = 2*ell+1
        qext[idx] = 2/x**2*np.sum(weight*np.real(a+b))
        qsca[idx] = 2/x**2*np.sum(weight*(np.abs(a)**2+np.abs(b)**2))
    qabs = qext-qsca
    if not np.all(np.isfinite(qext+qsca)) or np.any(qabs < -1e-8):
        raise ValueError("Mie numerical instability or nonphysical absorption")
    qabs = np.maximum(qabs, 0)
    result = {"Qext":qext, "Qsca":qsca, "Qabs":qabs,
              "Cext_nm2":qext*np.pi*radius**2,
              "Csca_nm2":qsca*np.pi*radius**2,
              "Cabs_nm2":qabs*np.pi*radius**2}
    return {k: v.item() if v.ndim == 0 else v for k,v in result.items()}
