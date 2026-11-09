import numpy as np
from nanoparticle_optics.film_tmm import film_rt
from nanoparticle_optics.preview_server import app

def test_lossless_film_energy():
    wl=np.linspace(450,900,21)
    out=film_rt(wl,np.full_like(wl,1.7+0j),250,1.5)
    np.testing.assert_allclose(out["R"]+out["T"],1,atol=1e-12)
    np.testing.assert_allclose(out["A"],0,atol=1e-12)

def test_absorbing_film_energy():
    wl=np.linspace(450,900,21)
    out=film_rt(wl,np.full_like(wl,1.7+.1j),250,1.5)
    assert np.all(out["A"]>=0)
    np.testing.assert_allclose(out["R"]+out["T"]+out["A"],1,atol=1e-12)

def test_zero_thickness_matches_bare_interface():
    wl=np.array([500.,700.])
    out=film_rt(wl,1.7+.1j,0,1.5)
    np.testing.assert_allclose(out["R"],((1-1.5)/(1+1.5))**2,atol=1e-12)

def test_api_returns_film_rta():
    response=app.test_client().post("/api/calculate",json={"film_thickness_nm":200,"substrate_n":1.5,"points":5})
    assert response.status_code==200,response.get_json()
    data=response.get_json()
    np.testing.assert_allclose(np.array(data["R"])+data["T"]+data["A"],1,atol=1e-9)
