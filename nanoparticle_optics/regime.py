"""Geometry-aware, conservative particle-scattering applicability diagnostics.

Diagnostics are advisory; they never silently substitute Rayleigh for Mie.
"""
import numpy as np

def scattering_regime(wavelength_nm,radius_nm,host_n,particle_n,distribution_radii_nm=None,aggregate=False,shape="sphere"):
    wl=np.asarray(wavelength_nm,float)
    nh=np.asarray(host_n,complex)
    npart=np.asarray(particle_n,complex)
    if wl.ndim!=1 or not len(wl) or np.any(~np.isfinite(wl)) or np.any(wl<=0):
        raise ValueError("Positive finite wavelength array required")
    if radius_nm<=0 or not np.isfinite(radius_nm):raise ValueError("Positive particle radius required")
    if np.any(nh.real<=0) or np.any(np.abs(nh.imag)>1e-12):
        raise ValueError("Regime guidance requires a lossless positive-index host")
    if shape not in ("sphere","spheroid","cylinder","irregular"):
        raise ValueError("Unknown particle shape")
    radii=np.asarray(distribution_radii_nm if distribution_radii_nm is not None else [radius_nm],float)
    if radii.ndim!=1 or radii.size==0 or np.any(~np.isfinite(radii)) or np.any(radii<=0):
        raise ValueError("Radii must be positive and finite")
    x=2*np.pi*nh.real[:,None]*radii[None,:]/wl[:,None]
    m=npart[:,None]/nh[:,None]
    mx=np.abs(m)*x
    flags=[]
    if shape!="sphere":
        flags.append("Homogeneous-sphere Mie outputs are not applicable to the selected nonspherical geometry; use an appropriate external solver.")
    if aggregate:
        flags.append("Aggregate coupling, structure factor and dependent scattering are not included in independent-sphere Mie.")
    if x.max()>=1:
        flags.append("Geometric or higher-order multipoles may contribute; do not use Rayleigh approximation across the entire sweep.")
    if np.max(mx)>=.3:
        flags.append("Rayleigh small-internal-phase criterion |m|x << 1 is not assured.")
    if not flags:flags.append("Small size parameters suggest a Rayleigh-like regime; full Mie remains the calculated model.")
    return {"x_min":float(x.min()),"x_max":float(x.max()),"mx_max":float(mx.max()),
            "rayleigh_heuristic_all":bool(np.max(x)<.3 and np.max(mx)<.3),
            "shape":shape,"aggregate":bool(aggregate),"warnings":flags,
            "solver":"independent_homogeneous_sphere_Mie" if shape=="sphere" and not aggregate else "external_solver_required"}
