import numpy as np
import pytest
from nanoparticle_optics.three_phase import three_phase_epsilon
from nanoparticle_optics.dispersion_models import model_epsilon
from nanoparticle_optics.preview_server import app

def test_three_phase_limits():
    h=np.array([2.25+0j,2.3+0j]);p=np.array([3.5+0j,3.6+0j]);v=np.ones(2)
    np.testing.assert_allclose(three_phase_epsilon(h,p,v,0,0),h)
    np.testing.assert_allclose(three_phase_epsilon(h,p,v,.1,0),
        three_phase_epsilon(h,p,v,.1,0,method="maxwell_garnett"))
    assert np.all(np.isfinite(three_phase_epsilon(h,p,v,.1,.2,method="bruggeman")))
    with pytest.raises(ValueError,match="sum"):
        three_phase_epsilon(h,p,v,.7,.4)

def test_complex_lorentz_passive():
    wl=np.linspace(450,800,15)
    for kind in ("lorentz","drude_lorentz"):
        model={"model":kind,"epsilon_inf":2.5,"oscillators":[{"strength":.8,"energy_ev":2.3,"gamma_ev":.12}]}
        if kind=="drude_lorentz":model.update(plasma_ev=1,drude_gamma_ev=.1)
        eps=model_epsilon(model,wl)
        assert np.all(eps.imag>=0)
        assert np.ptp(eps.real)>.01

def test_three_phase_api_energy():
    r=app.test_client().post("/api/calculate",json={"start_nm":500,"end_nm":800,
        "points":7,"fraction":.1,"pore_fraction":.2,"pore_n":1.33,
        "particle_model":{"model":"silica_malitson"},
        "host_model":{"model":"cauchy","A":1.48,"B_um2":.002},
        "film_thickness_nm":200,"substrate_n":1.5})
    assert r.status_code==200,r.get_json()
    d=r.get_json()
    assert abs(d["host_fraction"]-.7)<1e-12
    np.testing.assert_allclose(np.array(d["R"])+d["T"]+d["A"],1,atol=1e-9)
