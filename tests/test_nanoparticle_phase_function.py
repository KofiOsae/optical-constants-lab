import numpy as np
from nanoparticle_optics.phase_function import sphere_phase_function
from nanoparticle_optics.angular import mie_anisotropy

def test_phase_normalization_and_g():
    wl=550.;radius=65.;index=1.7+.02j
    p=sphere_phase_function(wl,radius,index,1.33,angles_deg=np.linspace(0,180,721))
    theta=np.deg2rad(p["angle_deg"])
    integral=2*np.pi*np.trapz(p["phase_sr_inv"]*np.sin(theta),theta)
    np.testing.assert_allclose(integral,1,rtol=1e-8)
    g=mie_anisotropy(np.array([wl]),radius,index,1.33)["g"][0]
    np.testing.assert_allclose(p["g_angular"],g,atol=5e-4)
    assert np.all(p["perpendicular_nm2_sr"]>=0)
    assert np.all(p["parallel_nm2_sr"]>=0)

def test_rayleigh_angular_symmetry():
    p=sphere_phase_function(600,3,1.5,1.,angles_deg=np.linspace(0,180,721))
    assert abs(p["g_angular"])<.05
