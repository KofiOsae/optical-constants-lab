"""Literature-grounded forward-model checks.

Malitson (1965), JOSA 55, 1205-1209: fused-silica Sellmeier coefficients,
reproduced in refractiveindex.info, database entry SiO2/Malitson.
These tests validate material dispersion and passive coherent propagation,
NOT blind inversion accuracy.
"""
import numpy as np
from optical_inversion_package.inversion import coherent_film

def malitson_silica(wavelength_um):
    l=np.asarray(wavelength_um,dtype=float)
    l2=l*l
    return np.sqrt(1+0.6961663*l2/(l2-0.0684043**2)
                     +0.4079426*l2/(l2-0.1162414**2)
                     +0.8974794*l2/(l2-9.896161**2))

def test_malitson_known_refractive_indices():
    # Published silica dispersion yields ~1.458 at 633 nm and ~1.444 at 1550 nm.
    assert abs(float(malitson_silica(.6328))-1.4570)<.002
    assert abs(float(malitson_silica(1.55))-1.4440)<.002

def test_literature_silica_lossless_film_energy_conservation():
    w=np.linspace(350,1400,106)
    n=malitson_silica(w/1000.)
    ns=malitson_silica(w/1000.)
    R,T=coherent_film(w,n.astype(complex),2000.,ns)
    assert np.all(np.isfinite(R)) and np.all(np.isfinite(T))
    assert np.max(np.abs(R+T-1))<1e-10
    # Identical film and substrate: only air/silica interface is optically present.
    expected=((1-n)/(1+n))**2
    np.testing.assert_allclose(R,expected,rtol=1e-10,atol=1e-12)
