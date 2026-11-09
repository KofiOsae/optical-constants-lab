import numpy as np
from nanoparticle_optics.preview_server import app

def test_effective_complex_permittivity_and_nk_consistency():
    client=app.test_client()
    r=client.post("/api/calculate",json={"host_n":1.5,"host_k":0,
      "particle_n":2.,"particle_k":0.1,"fraction":0.1,
      "radius_nm":30,"start_nm":500,"end_nm":700,"points":5})
    assert r.status_code==200,r.get_json()
    d=r.get_json()
    nk=np.array(d["effective_n"])+1j*np.array(d["effective_k"])
    eps=np.array(d["effective_epsilon1"])+1j*np.array(d["effective_epsilon2"])
    np.testing.assert_allclose(nk**2,eps,rtol=1e-11,atol=1e-11)
    assert np.all(eps.imag>=0)

def test_host_limit():
    r=app.test_client().post("/api/calculate",json={"fraction":0,"host_n":1.5,
      "host_k":0,"particle_n":2.,"particle_k":0.1,"points":4})
    assert r.status_code==200,r.get_json()
    np.testing.assert_allclose(r.get_json()["effective_epsilon1"],2.25,rtol=1e-11)
    np.testing.assert_allclose(r.get_json()["effective_epsilon2"],0,atol=1e-11)
