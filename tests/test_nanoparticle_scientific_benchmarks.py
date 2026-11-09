"""Analytic and convergence checks for isolated nanoparticle modules.

No production imports. Independent reference limits are derived from
electrostatic polarizability and the optical theorem.
"""
import numpy as np
import pytest
from nanoparticle_optics.mie import mie_sphere
from nanoparticle_optics.effective_medium import maxwell_garnett, bruggeman
from nanoparticle_optics.ensemble import ensemble_cross_sections, lognormal_radius_bins


def test_mie_rayleigh_absorbing_sphere():
    wl, radius, n = 600., 2., 1.5+0.05j
    x = 2*np.pi*radius/wl
    alpha = (n*n-1)/(n*n+2)
    qabs_limit = 4*x*alpha.imag
    actual = mie_sphere(wl, radius, n)
    assert np.isclose(actual["Qabs"], qabs_limit, rtol=.025)


def test_mie_rayleigh_scattering_size_scaling():
    wl, n = 600., 1.5
    q1 = mie_sphere(wl, 3., n)["Qsca"]
    q2 = mie_sphere(wl, 6., n)["Qsca"]
    assert np.isclose(q2/q1, 16., rtol=.02)


def test_mie_extinction_balance_and_nonnegative_absorption_grid():
    for radius in [10., 30., 100., 200.]:
        for index in [1.2+0j, 1.8+0.05j, 2.5+0.3j]:
            result = mie_sphere(np.array([400., 550., 800.]), radius, index)
            assert np.all(result["Qabs"] >= 0)
            assert np.allclose(result["Qext"], result["Qsca"]+result["Qabs"], atol=1e-8)


def test_mg_first_order_polarizability():
    host, inclusion = 2.25+0j, 4.+0.1j
    f = 1e-5
    predicted = host + 3*f*host*(inclusion-host)/(inclusion+2*host)
    assert np.isclose(maxwell_garnett(host, inclusion, f), predicted, rtol=1e-8)


def test_bruggeman_complex_equation_grid():
    h, i = 2.25+0.01j, 4.+0.2j
    for f in [.05, .2, .4, .6, .8, .95]:
        e = bruggeman(h, i, f)
        residual = f*(i-e)/(i+2*e)+(1-f)*(h-e)/(h+2*e)
        assert abs(residual) < 1e-10
        assert e.imag >= -1e-10


def test_lognormal_ensemble_convergence():
    wl = np.array([550.])
    outputs = []
    for bins in [81, 161, 321]:
        radii, weights = lognormal_radius_bins(30., 1.25, count=bins)
        ensemble = ensemble_cross_sections(wl, radii, weights, 1.65+0.03j)
        outputs.append(ensemble["Csca_nm2"][0])
    assert abs(outputs[-2]/outputs[-1]-1) < .001
    assert abs(outputs[0]/outputs[-1]-1) < .005


def test_ensemble_extinction_additivity():
    radii, weights = lognormal_radius_bins(40., 1.3, count=41)
    result = ensemble_cross_sections(np.array([450., 600.]), radii, weights, 1.8+0.1j)
    assert np.allclose(result["Cext_nm2"],
                       result["Csca_nm2"]+result["Cabs_nm2"], atol=1e-8)
