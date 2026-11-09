"""Research-only ellipsometric inversion using existing stack_models forward solver.

Fit a single homogeneous isotropic film with constant n,k across the selected
spectral window. Delta residual is wrapped by fitting complex rho instead of
subtracting angles. Multi-angle data accepted. No production route imported.
"""
import numpy as np
from scipy.optimize import least_squares
from stack_models import ellipsometry

def fit_ellipsometry(wavelength_nm, psi_deg, delta_deg, angles_deg,
                     substrate_n=1.5, initial=(1.7,0.02,100.),
                     lower=(1.01,0.,1.), upper=(5.,3.,3000.),
                     max_nfev=500):
    wl=np.asarray(wavelength_nm,float)
    angles=np.atleast_1d(np.asarray(angles_deg,float))
    psi=np.asarray(psi_deg,float);delta=np.asarray(delta_deg,float)
    if wl.ndim!=1 or len(wl)<4 or not np.all(np.isfinite(wl)) or np.any(wl<=0):
        raise ValueError("At least four positive finite wavelengths required")
    if np.any(~np.isfinite(angles)) or np.any((angles<=0)|(angles>=89)):
        raise ValueError("Angles must be between 0 and 89 degrees")
    if psi.shape != (len(angles),len(wl)) or delta.shape!=psi.shape:
        raise ValueError("Psi/Delta must be shaped (angles, wavelengths)")
    if not np.all(np.isfinite(psi+delta)) or np.any((psi<0)|(psi>90)):
        raise ValueError("Invalid ellipsometric observations")
    lo=np.asarray(lower,float);hi=np.asarray(upper,float);x0=np.asarray(initial,float)
    if lo.shape!=(3,) or hi.shape!=(3,) or x0.shape!=(3,) or np.any(lo>=hi) or np.any(x0<lo) or np.any(x0>hi):
        raise ValueError("Invalid n,k,thickness bounds or initial values")
    def forward(x):
        rows=[]
        for angle in angles:
            out=ellipsometry(wl,[{"n":float(x[0]),"k":float(x[1]),"thickness_nm":float(x[2])}],
                             substrate_n=substrate_n,angle_deg=float(angle))
            rows.append(np.tan(np.deg2rad(out["psi_deg"]))*np.exp(1j*np.deg2rad(out["delta_deg"])))
        return np.asarray(rows)
    observed=np.tan(np.deg2rad(psi))*np.exp(1j*np.deg2rad(delta))
    def residual(x):
        diff=forward(x)-observed
        return np.r_[diff.real.ravel(),diff.imag.ravel()]
    fit=least_squares(residual,x0,bounds=(lo,hi),max_nfev=max_nfev)
    predicted=forward(fit.x)
    return {"n":float(fit.x[0]),"k":float(fit.x[1]),
            "thickness_nm":float(fit.x[2]),"success":bool(fit.success),
            "message":fit.message,"cost":float(fit.cost),
            "rmse_rho":float(np.sqrt(np.mean(abs(predicted-observed)**2))),
            "jacobian_rank":int(np.linalg.matrix_rank(fit.jac)),
            "nfev":int(fit.nfev)}
