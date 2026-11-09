import numpy as np
from nanoparticle_optics.angular import mie_anisotropy,ensemble_anisotropy
from nanoparticle_optics.mie import mie_sphere

def test_rayleigh_g_tends_to_zero():
    v=mie_anisotropy(np.array([600.]),3,1.5,1.)
    assert abs(v["g"][0])<0.05

def test_single_sphere_matches_mie_qsca():
    wl=np.array([450.,650.])
    v=mie_anisotropy(wl,45,1.7+.03j,1.33)
    m=mie_sphere(wl,45,1.7+.03j,1.33)
    np.testing.assert_allclose(v["Qsca"],m["Qsca"],rtol=1e-11)

def test_ensemble_transport_scattering_weighted():
    wl=np.array([500.,650.])
    radii=np.array([20.,80.])
    weights=np.array([.8,.2])
    density=1e-8
    e=ensemble_anisotropy(wl,radii,weights,1.7+.02j,1.33,density)
    parts=[mie_anisotropy(wl,r,1.7+.02j,1.33) for r in radii]
    c=sum(w*p["Csca_nm2"] for w,p in zip(weights,parts))
    moment=sum(w*p["Csca_nm2"]*p["g"] for w,p in zip(weights,parts))
    np.testing.assert_allclose(e["g"],moment/c)
    np.testing.assert_allclose(e["mu_s_prime_nm_inv"],density*(c-moment))
    assert np.all(e["g"]<=1) and np.all(e["g"]>=-1)
