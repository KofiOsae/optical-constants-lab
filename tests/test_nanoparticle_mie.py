import numpy as np
import pytest
from nanoparticle_optics.mie import mie_sphere

def test_identical_indices_have_zero_scattering():
    out = mie_sphere([400,600,800], 60, 1.5, 1.5)
    assert np.max(np.abs(out["Qsca"])) < 1e-12
    assert np.max(np.abs(out["Qext"])) < 1e-12

def test_lossless_energy_conservation():
    out = mie_sphere([400,500,800], 80, 1.8, 1.0)
    assert np.allclose(out["Qext"], out["Qsca"], rtol=1e-9, atol=1e-11)
    assert np.all(out["Qabs"] >= 0)

def test_absorbing_sphere_passive():
    out = mie_sphere([450,550,650], 40, 2+0.2j, 1.0)
    assert np.all(out["Qabs"] > 0)
    assert np.allclose(out["Qext"], out["Qsca"]+out["Qabs"], atol=1e-10)

def test_rayleigh_limit():
    wavelength, radius, n = 600., 3., 1.5
    x = 2*np.pi*radius/wavelength
    m = n
    rayleigh = (8/3)*x**4*abs((m*m-1)/(m*m+2))**2
    actual = mie_sphere(wavelength, radius, n)["Qsca"]
    assert np.isclose(actual, rayleigh, rtol=.01)

def test_cross_section_consistency():
    out = mie_sphere(500, 40, 1.7)
    assert np.isclose(out["Csca_nm2"], out["Qsca"]*np.pi*40**2)

def test_invalid_inputs():
    with pytest.raises(ValueError):
        mie_sphere(0, 10, 1.5)
    with pytest.raises(ValueError):
        mie_sphere(500, 10, 1.5-0.1j)
