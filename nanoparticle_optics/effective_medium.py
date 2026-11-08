"""Independent quasistatic effective-medium models (not wired to production).

Time convention: exp(-i omega t); passive media have Im(epsilon) >= 0.
Volume fraction refers to inclusions. Models assume isotropic, spherical,
subwavelength inclusions and neglect interparticle multiple scattering.
"""
import numpy as np


def _inputs(eps_host, eps_inclusion, fraction):
    host, inc, f = np.broadcast_arrays(
        np.asarray(eps_host, dtype=complex),
        np.asarray(eps_inclusion, dtype=complex),
        np.asarray(fraction, dtype=float),
    )
    if not np.all(np.isfinite(host)) or not np.all(np.isfinite(inc)):
        raise ValueError("Permittivities must be finite")
    if not np.all(np.isfinite(f)) or np.any((f < 0) | (f > 1)):
        raise ValueError("Volume fraction must be in [0, 1]")
    return host, inc, f


def maxwell_garnett(eps_host, eps_inclusion, fraction):
    """Spherical inclusions in a connected host, quasistatic approximation."""
    h, i, f = _inputs(eps_host, eps_inclusion, fraction)
    numerator = i + 2*h + 2*f*(i-h)
    denominator = i + 2*h - f*(i-h)
    with np.errstate(divide="ignore", invalid="ignore"):
        result = h * numerator / denominator
    result = np.where(f == 0, h, np.where(f == 1, i, result))
    if not np.all(np.isfinite(result)):
        raise ValueError("Maxwell-Garnett singularity: check material parameters")
    return result.item() if result.ndim == 0 else result


def bruggeman(eps_host, eps_inclusion, fraction):
    """Symmetric two-component spherical Bruggeman EMA.

    Solves 2*e**2 - B*e - h*i = 0 with
    B=(2-3*f)*h+(3*f-1)*i. Select the root closest to the
    dilute Maxwell-Garnett prediction; ambiguous near resonances
    should be checked for continuity and passivity.
    """
    h, i, f = _inputs(eps_host, eps_inclusion, fraction)
    b = (2-3*f)*h + (3*f-1)*i
    disc = np.sqrt(b*b + 8*h*i)
    r1, r2 = (b+disc)/4, (b-disc)/4
    reference = maxwell_garnett(h, i, f)
    result = np.where(np.abs(r1-reference) <= np.abs(r2-reference), r1, r2)
    result = np.where(f == 0, h, np.where(f == 1, i, result))
    if not np.all(np.isfinite(result)):
        raise ValueError("Bruggeman singularity: check material parameters")
    return result.item() if result.ndim == 0 else result
