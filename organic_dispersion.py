"""Passive, causal oscillator lineshapes for organic molecular films.

Photon energy E is in eV, wavelength in nm; dielectric convention exp(-i omega t).
All susceptibilities are positive sums of damped Lorentz poles. This ensures a
nonnegative imaginary dielectric part for positive photon energies and
Kramers-Kronig compatibility, unlike directly fitting arbitrary k Gaussians.
"""
from functools import lru_cache
from math import exp, factorial
import numpy as np
from scipy.special import wofz

HC_EV_NM = 1239.841984

MODEL_PARAMETER_COUNTS = {
    "modified_lorentz": 4,  # f, E0, width, broad_component_fraction
    "kk_gaussian": 3,       # f, E0, sigma
    "voigt": 4,             # f, E0, gamma, sigma
    "vibronic": 5,          # f, E00, vibration_energy, gamma, Huang_Rhys
}

@lru_cache(maxsize=4)
def _gauss_hermite(order):
    nodes, weights = np.polynomial.hermite.hermgauss(order)
    return nodes, weights / np.sqrt(np.pi)

def _lorentz(E, strength, center, gamma):
    if strength < 0 or center <= 0 or gamma <= 0:
        raise ValueError("Oscillator strengths must be nonnegative, energies and widths positive")
    return strength * center**2 / (center**2 - E**2 - 1j * gamma * E)

def _inhomogeneous(E, strength, center, gamma, sigma, order=24):
    """Gaussian distribution of positive resonance energies; each member causal."""
    if sigma < 0:
        raise ValueError("Gaussian broadening sigma must be nonnegative")
    if sigma == 0:
        return _lorentz(E, strength, center, gamma)
    nodes, weights = _gauss_hermite(order)
    centers = center + np.sqrt(2) * sigma * nodes
    keep = centers > 0
    if not np.any(keep):
        raise ValueError("Gaussian distribution has no positive resonance energies")
    weights = weights[keep] / np.sum(weights[keep])
    out = np.zeros_like(E, dtype=complex)
    for c, w in zip(centers[keep], weights):
        out += w * _lorentz(E, strength, float(c), gamma)
    return out

def _analytic_gaussian_lorentz(E,strength,center,gamma,sigma):
    """Exact Gaussian convolution over signed resonance energy via Faddeeva.

    The squared resonance denominator makes negative centers passive; this
    symmetric-energy Gaussian extension is numerically stable for broad bands.
    Unlike fixed-node quadrature, it cannot create discretization sub-peaks.
    """
    if strength<0 or center<=0 or gamma<=0 or sigma<0:
        raise ValueError('Oscillator requires passive, positive energy/width parameters')
    if sigma==0:
        return _lorentz(E,strength,center,gamma)
    z=np.sqrt(np.asarray(E,float)**2+1j*gamma*np.asarray(E,float))
    spread=np.sqrt(2.)*sigma
    minus=1j*np.sqrt(np.pi)/spread*wofz((z-center)/spread)
    plus=-1j*np.sqrt(np.pi)/spread*wofz((z+center)/spread)
    return strength*(1.+.5*z*(minus-plus))


def _tauc_lorentz(E, A, Eg, separation, C):
    """Jellison-Modine imaginary dielectric function and finite quadrature KK real part.

    E0=Eg+separation ensures the oscillator peak lies above the optical onset.
    KK integral uses 512-point Gauss-Legendre quadrature over [Eg, 80 eV];
    the finite tail truncation is an approximation, not an exact analytic TL.
    """
    from numpy.polynomial.legendre import leggauss
    if not (A >= 0 and Eg > 0 and separation > 0 and C > 0):
        raise ValueError('Tauc-Lorentz requires A>=0, Eg>0, E0>Eg and C>0')
    E=np.asarray(E,float)
    E0=Eg+separation
    def im(x):
        x=np.asarray(x,float)
        out=np.zeros_like(x)
        valid=x>Eg
        z=x[valid]
        out[valid]=A*C*E0*(z-Eg)**2/(z*((z*z-E0*E0)**2+C*C*z*z))
        return out
    x0,weights=_tl_nodes()
    upper=max(80., E0+30*C, float(np.max(E))*4.)
    x=Eg+(upper-Eg)*(x0+1)/2
    dx=weights*(upper-Eg)/2
    f=im(x)
    fp=im(E)
    denominator=x[None,:]**2-E[:,None]**2
    term=(f[None,:]-fp[:,None])*x[None,:]/np.where(np.abs(denominator)<1e-10,1e-10,denominator)
    integral=np.sum(term*dx[None,:],axis=1)
    # Analytic principal value for the subtracted constant eps2(E).
    lower_term=np.maximum(np.abs(Eg*Eg-E*E),1e-30)
    upper_term=np.maximum(np.abs(upper*upper-E*E),1e-30)
    integral+=.5*fp*np.log(upper_term/lower_term)
    return (2/np.pi)*integral+1j*fp

@lru_cache(maxsize=1)
def _tl_nodes():
    return np.polynomial.legendre.leggauss(512)

def organic_dielectric(wavelength_nm, parameters, model, oscillators=1, quadrature_order=24):
    """Evaluate dielectric function for a selected organic model.

    Parameter layout is [epsilon_inf, oscillator_1..., oscillator_2..., ...].
    Modified Lorentz: a positive mixture of a narrow and a broad Lorentz pole,
    both centered on E0. KK Gaussian is the narrow homogeneous-width limit
    of Gaussian distributed Lorentz oscillators. Voigt includes both widths.
    Vibronic uses nonnegative Poisson/Huang-Rhys weights for five sidebands.
    These are explicitly defined phenomenological model families, not claims
    of a unique microscopic molecular Hamiltonian.
    """
    w = np.asarray(wavelength_nm, dtype=float)
    p = np.asarray(parameters, dtype=float)
    if model not in MODEL_PARAMETER_COUNTS:
        raise ValueError("Unsupported organic dispersion model")
    if not 1 <= oscillators <= 8:
        raise ValueError("Oscillator count must be 1 through 8")
    expected = 1 + oscillators * MODEL_PARAMETER_COUNTS[model]
    if p.shape != (expected,) or not np.all(np.isfinite(p)):
        raise ValueError(f"{model} requires {expected} finite coefficients")
    if np.any(~np.isfinite(w)) or np.any(w <= 0) or p[0] <= 0:
        raise ValueError("Wavelengths and epsilon_inf must be finite and positive")
    E = HC_EV_NM / w
    eps = np.full(w.shape, p[0], dtype=complex)
    stride = MODEL_PARAMETER_COUNTS[model]
    for j in range(oscillators):
        q = p[1+j*stride:1+(j+1)*stride]
        if model == "modified_lorentz":
            f, e0, width, fraction = q
            if not 0 <= fraction <= 1:
                raise ValueError("Broad-component fraction must be between zero and one")
            eps += (1-fraction)*_lorentz(E,f,e0,width)
            eps += fraction*_lorentz(E,f,e0,3*width)
        elif model == "kk_gaussian":
            f,e0,sigma = q
            eps += _analytic_gaussian_lorentz(E,f,e0,0.005,sigma)
        elif model == "voigt":
            f,e0,gamma,sigma = q
            eps += _inhomogeneous(E,f,e0,gamma,sigma,quadrature_order)
        elif model == "tauc_lorentz":
            eps += _tauc_lorentz(E,*q)
        else:
            f,e00,delta,gamma,S = q
            if delta <= 0 or S < 0:
                raise ValueError("Vibrational spacing must be positive and Huang-Rhys factor nonnegative")
            weights = np.array([exp(-S)*S**m/factorial(m) for m in range(6)])
            weights /= weights.sum()
            for m,weight in enumerate(weights):
                eps += _lorentz(E,f*weight,e00+m*delta,gamma)
    return eps

def organic_nk(wavelength_nm, parameters, model, oscillators=1):
    eps = organic_dielectric(wavelength_nm, parameters, model, oscillators)
    N = np.sqrt(eps + 0j)
    return N.real, N.imag

def organic_parameter_setup(model, oscillators, low_wavelength_nm, high_wavelength_nm):
    """Defaults and physically admissible bounds shared by all inverse fitters."""
    if model not in MODEL_PARAMETER_COUNTS or not 1 <= oscillators <= 8:
        raise ValueError("Unsupported organic model or oscillator count")
    elo=HC_EV_NM/high_wavelength_nm
    ehi=HC_EV_NM/low_wavelength_nm
    starts=[2.5]; lower=[1.]; upper=[25.]
    for energy in np.linspace(elo,ehi,oscillators+2)[1:-1]:
        if model=="modified_lorentz":
            values=[.3,float(energy),.20,.35]
            lb=[0.,.05,.001,0.]
            ub=[20.,max(10.,2*ehi),5.,1.]
        elif model=="kk_gaussian":
            values=[.3,float(energy),.15]
            lb=[0.,.05,.001]
            ub=[20.,max(10.,2*ehi),5.]
        elif model=="voigt":
            values=[.3,float(energy),.08,.12]
            lb=[0.,.05,.001,0.]
            ub=[20.,max(10.,2*ehi),5.,5.]
        elif model=="tauc_lorentz":
            values=[3.,max(.2,float(energy)*.65),max(.1,float(energy)*.35),.3]
            lb=[0.,.02,.02,.005]
            ub=[100.,max(10.,2*ehi),max(10.,2*ehi),8.]
        else:
            values=[.3,float(energy),.17,.12,1.]
            lb=[0.,.05,.001,.001,0.]
            ub=[20.,max(10.,2*ehi),2.,5.,5.]
        starts.extend(values);lower.extend(lb);upper.extend(ub)
    return np.array(starts,float),np.array(lower,float),np.array(upper,float)
