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
