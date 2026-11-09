"""Independent particle-ensemble optics, dilute independent-scattering limit.

Inputs are discrete radius bins with *number* probability weights, not
volume probabilities. All cross sections are in nm^2; number density in
particles/nm^3. Extinction coefficient is in inverse nm. This model does
not include dependent scattering, aggregation, or thin-film interference.
"""
import numpy as np
from nanoparticle_optics.mie import mie_sphere


def lognormal_radius_bins(median_radius_nm, geometric_std, count=81, sigmas=4):
    """Return bin centers and normalized number probabilities.

    Logarithmic bin boundaries are used so that the weights represent
    probability mass rather than samples of the probability density.
    """
    if not np.isfinite(median_radius_nm) or median_radius_nm <= 0:
        raise ValueError("Median radius must be positive")
    if not np.isfinite(geometric_std) or geometric_std < 1:
        raise ValueError("Geometric standard deviation must be >= 1")
    if not isinstance(count, int) or count < 3:
        raise ValueError("At least three bins required")
    if not np.isfinite(sigmas) or sigmas <= 0:
        raise ValueError("sigmas must be positive")
    if geometric_std == 1:
        return np.array([float(median_radius_nm)]), np.array([1.])
    from scipy.special import ndtr
    s = np.log(geometric_std)
    edges = np.linspace(-sigmas, sigmas, count+1)
    radii = median_radius_nm * np.exp(s*(edges[:-1]+edges[1:])/2)
    probabilities = np.diff(ndtr(edges))
    probabilities /= probabilities.sum()
    return radii, probabilities


def ensemble_cross_sections(wavelength_nm, radii_nm, number_weights,
                            n_particle, n_host=1.0,
                            number_density_nm3=None):
    """Average independent-sphere cross sections over number probabilities.

    Returns number-averaged Cext/Csca/Cabs, mean particle volume and,
    optionally, macroscopic attenuation coefficients (1/nm).
    n_particle may be a scalar or a 1D wavelength-indexed array.
    """
    wl = np.asarray(wavelength_nm, dtype=float)
    r = np.asarray(radii_nm, dtype=float)
    weights = np.asarray(number_weights, dtype=float)
    if wl.ndim != 1 or not len(wl) or np.any(~np.isfinite(wl)) or np.any(wl <= 0):
        raise ValueError("Wavelengths must be a positive finite 1D array")
    if r.ndim != 1 or weights.shape != r.shape or not len(r):
        raise ValueError("Radius and weight arrays must be matching 1D arrays")
    if np.any(~np.isfinite(r)) or np.any(r <= 0) or np.any(~np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("Invalid radii or number weights")
    total = weights.sum()
    if total <= 0:
        raise ValueError("Number weights must have positive sum")
    weights = weights/total
    index = np.asarray(n_particle, dtype=complex)
    if index.ndim > 1 or (index.ndim == 1 and index.shape != wl.shape):
        raise ValueError("Particle index must be scalar or wavelength-indexed")
    if number_density_nm3 is not None:
        density = float(number_density_nm3)
        if not np.isfinite(density) or density < 0:
            raise ValueError("Number density must be nonnegative")
    else:
        density = None
    mean = {k: np.zeros(wl.shape) for k in ("Cext_nm2", "Csca_nm2", "Cabs_nm2")}
    for radius, weight in zip(r, weights):
        if weight == 0:
            continue
        result = mie_sphere(wl, radius, index, n_host)
        for key in mean:
            mean[key] += weight*result[key]
    mean_volume = float(np.sum(weights*(4*np.pi/3)*r**3))
    result = {"wavelength_nm":wl, **mean, "mean_particle_volume_nm3":mean_volume,
              "number_weights":weights, "volume_weights":weights*r**3/np.sum(weights*r**3)}
    if density is not None:
        result["volume_fraction"] = density*mean_volume
        if result["volume_fraction"] > 1:
            raise ValueError("Particle volume fraction cannot exceed 1")
        for key in mean:
            result[key.replace("C", "mu_").replace("_nm2", "_nm_inv")] = density*mean[key]
    return result
