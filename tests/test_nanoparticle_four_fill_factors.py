import numpy as np
from nanoparticle_optics.preview_server import app

def test_four_fill_factor_ema_silica_pmma():
    client=app.test_client()
    results=[]
    for f in (0,.05,.10,.20):
        r=client.post("/api/calculate",json={
            "method":"maxwell_garnett","fraction":f,"start_nm":500,
            "end_nm":800,"points":7,"radius_nm":30,
            "host_model":{"model":"cauchy","A":1.48,"B_um2":.002},
            "particle_model":{"model":"silica_malitson"},
            "film_thickness_nm":200,"substrate_n":1.5})
        assert r.status_code==200,r.get_json()
        results.append(r.get_json())
    for d in results:
        eps=np.asarray(d["effective_epsilon1"])+1j*np.asarray(d["effective_epsilon2"])
        nk=np.asarray(d["effective_n"])+1j*np.asarray(d["effective_k"])
        np.testing.assert_allclose(nk**2,eps,atol=1e-10)
        np.testing.assert_allclose(np.asarray(d["R"])+d["T"]+d["A"],1,atol=1e-9)
    np.testing.assert_allclose(results[0]["effective_n"],results[0]["host_n"],atol=1e-12)
    assert all(len(d["effective_n"])==7 for d in results)
