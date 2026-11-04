"""Independent normal-incidence Fresnel benchmark for v2 teaching CSVs.

No imports of the Optical Constants Laboratory forward/inverse engines.
Checks generated datasets against an independent recursive amplitude formula.
This is *forward* verification, not blind inversion or experimental validation.
Run: python -m pytest -q tests/test_independent_fresnel_reference.py
"""
from pathlib import Path

import numpy as np
import pytest

DATA = Path(__file__).resolve().parents[1] / "examples" / "v2"


def independent_rt(wavelength_nm, layers, substrate_n=1.46, ambient_n=1.0):
    """Complex-amplitude recursion; layers listed incident-side first.

    r_j=(r_j,j+1+r_next*p^2)/(1+r_j,j+1*r_next*p^2)
    t_j=t_j,j+1*t_next*p/(1+r_j,j+1*r_next*p^2)
    p=exp(i*2*pi*N*d/lambda).
    R=|r|^2, T=Re(N_substrate)/Re(N_ambient)*|t|^2 for real
    positive ambient and substrate indices.
    """
    w = np.asarray(wavelength_nm, dtype=float)
    if w.ndim != 1 or not np.all(w > 0):
        raise ValueError("Wavelengths must be a positive one-dimensional array")
    media = [complex(ambient_n)] + [np.asarray(n, dtype=complex) for n, _ in layers] + [complex(substrate_n)]
    r_eff = np.zeros_like(w, dtype=complex)
    t_eff = np.ones_like(w, dtype=complex)
    for j in range(len(layers), -1, -1):
        n_left, n_right = media[j], media[j + 1]
        r = (n_left - n_right) / (n_left + n_right)
        t = 2 * n_left / (n_left + n_right)
        if j == len(layers):
            r_eff = np.broadcast_to(r, w.shape).astype(complex).copy()
            t_eff = np.broadcast_to(t, w.shape).astype(complex).copy()
        else:
            thickness = layers[j][1]
            phase = np.exp(2j * np.pi * n_right * thickness / w)
            denominator = 1 + r * r_eff * phase**2
            r_eff, t_eff = (
                (r + r_eff * phase**2) / denominator,
                (t * t_eff * phase) / denominator,
            )
    return np.abs(r_eff)**2, (substrate_n / ambient_n) * np.abs(t_eff)**2


@pytest.mark.parametrize("filename,geometry", [
    ("transparent_film_120nm.csv", "transparent120"),
    ("transparent_film_240nm.csv", "transparent240"),
    ("absorbing_film_180nm.csv", "absorbing180"),
    ("multilayer_cap50_film150.csv", "multilayer"),
])
def test_independent_fresnel_reference(filename, geometry):
    spectrum = np.genfromtxt(DATA / filename, delimiter=",", names=True)
    w = spectrum["wavelength_nm"]
    assert len(w) == 81
    assert np.all(np.diff(w) > 0)
    if geometry.startswith("transparent"):
        thickness = 120.0 if geometry == "transparent120" else 240.0
        layers = [(1.65, thickness)]
    elif geometry == "absorbing180":
        n = 1.73 + 0.12 * (550.0 / w - 1.0)
        k = 0.07 * np.exp(-0.5 * ((w - 610.0) / 60.0)**2)
        layers = [(n + 1j*k, 180.0)]
    else:
        layers = [(1.80, 50.0), (1.65, 150.0)]
    ref_R, ref_T = independent_rt(w, layers)
    err_R = float(np.max(np.abs(ref_R - spectrum["R"])))
    err_T = float(np.max(np.abs(ref_T - spectrum["T"])))
    assert err_R < 1e-8, (filename, "R", err_R)
    assert err_T < 1e-8, (filename, "T", err_T)
    assert np.all(spectrum["R"] >= 0)
    assert np.all(spectrum["T"] >= 0)
    assert np.max(spectrum["R"] + spectrum["T"]) <= 1 + 1e-8
    if geometry != "absorbing180":
        assert np.max(np.abs(spectrum["R"] + spectrum["T"] - 1)) < 1e-8


def test_absorbing_case_dissipates_energy():
    a = np.genfromtxt(DATA / "absorbing_film_180nm.csv", delimiter=",", names=True)
    absorption = 1.0 - a["R"] - a["T"]
    assert np.all(absorption >= -1e-8)
    assert float(np.max(absorption)) > 0.01
