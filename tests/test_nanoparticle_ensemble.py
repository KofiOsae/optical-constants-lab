import numpy as np
import pytest
from nanoparticle_optics.ensemble import lognormal_radius_bins, ensemble_cross_sections
from nanoparticle_optics.mie import mie_sphere

def test_monodisperse_equals_single_sphere():
    wl = np.array([450., 550., 650.])
    a = ensemble_cross_sections(wl, [40.], [1.], 1.7+0.05j)
    b = mie_sphere(wl, 40., 1.7+0.05j)
    for key in ("Cext_nm2", "Csca_nm2", "Cabs_nm2"):
        assert np.allclose(a[key], b[key], rtol=1e-12)

def test_weighted_average_and_volume_weights():
    wl = np.array([500.])
    r = [20., 40.]
    a = ensemble_cross_sections(wl, r, [1, 3], 1.6)
    c0 = mie_sphere(wl, 20., 1.6)["Csca_nm2"]
    c1 = mie_sphere(wl, 40., 1.6)["Csca_nm2"]
    assert np.allclose(a["Csca_nm2"], .25*c0+.75*c1)
    assert np.allclose(a["volume_weights"].sum(), 1)
    assert a["volume_weights"][1] > .75

def test_number_density_and_volume_fraction():
    wl = np.array([500.])
    a = ensemble_cross_sections(wl, [10.], [1.], 1.5, number_density_nm3=1e-7)
    assert np.allclose(a["mu_ext_nm_inv"], a["Cext_nm2"]*1e-7)
    assert np.isclose(a["volume_fraction"], 1e-7*4*np.pi/3*10**3)

def test_lognormal_probabilities():
    r, p = lognormal_radius_bins(30., 1.3, count=101)
    assert len(r) == len(p) == 101
    assert np.isclose(p.sum(), 1)
    assert np.all(p >= 0)
    assert np.all(np.diff(r) > 0)

def test_degenerate_distribution():
    r, p = lognormal_radius_bins(30., 1.)
    assert np.array_equal(r, [30.])
    assert np.array_equal(p, [1.])

def test_invalid_inputs():
    with pytest.raises(ValueError):
        ensemble_cross_sections([500.], [20.], [-1.], 1.5)
    with pytest.raises(ValueError):
        ensemble_cross_sections([500.], [20.], [1.], 1.5, number_density_nm3=1.)
