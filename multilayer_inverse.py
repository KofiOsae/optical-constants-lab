"""Constrained single-unknown-film inverse fit in a coherent known-layer stack.

Layer order: ambient -> before -> unknown -> after -> substrate.
The forward model is imported, not duplicated. Fixed known-layer indices
may be scalar or arrays sampled on the measurement wavelength grid.
"""
from dataclasses import dataclass
import numpy as np
from scipy.optimize import least_squares
from multilayer_forward import coherent_stack

@dataclass
class StackSpectrum:
    wavelength_nm: np.ndarray
    R: np.ndarray | None
    T: np.ndarray | None
    before: tuple = ()
    after: tuple = ()
    substrate_n: object = 1.5
    ambient_n: object = 1.
    angle_deg: float = 0.
    polarization: str = 's'
    sigma_R: float = .01
    sigma_T: float = .01

def predict_stack(sample, unknown_n, unknown_k, thickness_nm):
    layers=list(sample.before)+[(np.asarray(unknown_n)+1j*np.asarray(unknown_k),thickness_nm)]+list(sample.after)
    return coherent_stack(sample.wavelength_nm,layers,sample.substrate_n,
                          sample.ambient_n,sample.angle_deg,sample.polarization)

def fit_constant_unknown(sample, thickness_nm, n0=1.8, k0=.02,
                         fit_thickness=False, thickness_bounds=None,
                         max_nfev=300):
    """Proof-of-inversion benchmark for a constant-index unknown film.

    This is a bounded diagnostic baseline, not the final dispersive
    n(lambda),k(lambda) production fitting interface.
    """
    if sample.R is None and sample.T is None:
        raise ValueError('At least R or T is required')
    w=np.asarray(sample.wavelength_nm,float)
    if w.ndim!=1 or len(w)<3:raise ValueError('At least three wavelengths required')
    if thickness_nm<=0:raise ValueError('Thickness must be positive')
    if fit_thickness:
        if thickness_bounds is None or not (0<thickness_bounds[0]<thickness_bounds[1]):
            raise ValueError('Specify positive ordered thickness bounds')
        lo=[.5,0.,thickness_bounds[0]]
        hi=[6.,5.,thickness_bounds[1]]
        x0=[n0,k0,thickness_nm]
    else:
        lo=[.5,0.]
        hi=[6.,5.]
        x0=[n0,k0]
    def residual(x):
        r,t=predict_stack(sample,x[0],x[1],x[2] if fit_thickness else thickness_nm)
        parts=[]
        for measured,modeled,sigma in ((sample.R,r,sample.sigma_R),(sample.T,t,sample.sigma_T)):
            if measured is not None:
                data=np.asarray(measured,float)
                if data.shape!=w.shape or not np.all(np.isfinite(data)):
                    raise ValueError('Measurement must match finite wavelength grid')
                if sigma<=0:raise ValueError('Noise scale must be positive')
                parts.append((modeled-data)/sigma)
        return np.concatenate(parts)
    fit=least_squares(residual,x0,bounds=(lo,hi),max_nfev=max_nfev)
    return {'result':fit,'n':float(fit.x[0]),'k':float(fit.x[1]),
            'thickness_nm':float(fit.x[2] if fit_thickness else thickness_nm),
            'weighted_rmse':float(np.sqrt(np.mean(fit.fun**2)))}
