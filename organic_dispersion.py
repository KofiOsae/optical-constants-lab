"""Passive, causal oscillator lineshapes for organic molecular films.

Photon energy E is in eV, wavelength in nm; dielectric convention exp(-i omega t).
All susceptibilities are positive sums of damped Lorentz poles. This ensures a
nonnegative imaginary dielectric part for positive photon energies and
Kramers-Kronig compatibility, unlike directly fitting arbitrary k Gaussians.
"""
from functools import lru_cache
from math import exp, factorial
import numpy as np

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
            eps += _inhomogeneous(E,f,e0,0.005,sigma,quadrature_order)
        elif model == "voigt":
            f,e0,gamma,sigma = q
            eps += _inhomogeneous(E,f,e0,gamma,sigma,quadrature_order)
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
