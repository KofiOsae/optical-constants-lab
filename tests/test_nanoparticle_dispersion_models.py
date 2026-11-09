import numpy as np
import pytest
from nanoparticle_optics.dispersion_models import model_epsilon
from nanoparticle_optics.preview_server import app

def test_malitson_visible_index():
    eps=model_epsilon({"model":"silica_malitson"},[500,600,700])
    n=np.sqrt(eps)
    assert 1.45<n[1].real<1.47
    assert np.allclose(n.imag,0)
    assert n[0].real>n[-1].real

def test_cauchy_coefficients():
    eps=model_epsilon({"model":"cauchy","A":1.45,"B_um2":0.005},[500,700])
    assert eps[0].real>eps[-1].real

def test_silica_range_guard():
    with pytest.raises(ValueError,match="validity"):
        model_epsilon({"model":"silica_malitson"},[100,600])

def test_api_dispersion_model():
    r=app.test_client().post("/api/calculate",json={
        "particle_model":{"model":"silica_malitson"},
        "host_model":{"model":"cauchy","A":1.33,"B_um2":0.002},
        "start_nm":500,"end_nm":700,"points":3,"radius_nm":30,"fraction":0.05})
    assert r.status_code==200,r.get_json()
    d=r.get_json()
    assert d["particle_n"][0]>d["particle_n"][-1]
    assert len(d["effective_epsilon1"])==3
