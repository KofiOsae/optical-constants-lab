"""Coherent-film R/T import via the same validation used by the preview."""
import numpy as np
from optical_inversion_package.inversion import Spectrum
from material_import import inspect_material

def parse_spectrum(filename, content, thickness_nm, substrate_n=1., wavelength_unit='auto', signal_unit='auto'):
    data = inspect_material(filename, content, wavelength_unit, signal_unit)
    columns = data['columns']
    if 'R' not in columns and 'T' not in columns:
        raise ValueError('Inverse fitting requires reflectance R and/or transmittance T; this file can be previewed but not fitted with the current R/T model')
    w = np.asarray(data['wavelength_nm'], dtype=float)
    if len(w) < 8:
        raise ValueError('Inverse fitting requires at least 8 spectral points')
    R = np.asarray(columns['R'], dtype=float) if 'R' in columns else None
    T = np.asarray(columns['T'], dtype=float) if 'T' in columns else None
    d, ns = float(thickness_nm), float(substrate_n)
    if not np.isfinite(d) or not 1 <= d <= 1e7:
        raise ValueError('Film thickness must be between 1 and 10,000,000 nm')
    if not np.isfinite(ns) or not 1 <= ns <= 5:
        raise ValueError('Substrate refractive index must be between 1 and 5')
    return Spectrum(w, R, T, d, substrate_n=ns, name=str(filename or 'uploaded spectrum'))

def apply_rt_noise(spectrum, sigma=0., seed=2026, sample_index=0):
    """Add reproducible absolute Gaussian R/T noise ONCE before fitting.

    sigma is standard deviation in fractional R/T units, e.g. .002 means
    0.2 percentage points. Independent R/T realizations; seed+sample_index
    selects a stable per-sample stream. Input is copied to prevent mutation
    of imported/pristine spectra.
    """
    from dataclasses import replace
    sigma=float(sigma)
    seed=int(seed)
    if not np.isfinite(sigma) or sigma<0 or sigma>.1:
        raise ValueError('Noise standard deviation must be 0..0.1 in R/T fractions')
    if seed<0 or seed>2147483647 or sample_index<0 or sample_index>100:
        raise ValueError('Noise random seed or sample index outside supported range')
    rng=np.random.default_rng(np.random.SeedSequence([seed,int(sample_index)]))
    changes={}
    for ch in ('R','T'):
        observed=getattr(spectrum,ch)
        if observed is not None:
            values=np.array(observed,float,copy=True)
            if sigma:
                values=np.clip(values+rng.normal(0.,sigma,values.size),0.,1.)
            changes[ch]=values
    return replace(spectrum,**changes)
