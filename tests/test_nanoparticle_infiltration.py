import numpy as np
from nanoparticle_optics.preview_server import app

def test_pore_infiltration_changes_effective_index_and_film_rt():
    client=app.test_client()
    base={"method":"maxwell_garnett","fraction":.08,"pore_fraction":.2,
          "host_model":{"model":"cauchy","A":1.48,"B_um2":.002},
          "particle_model":{"model":"silica_malitson"},
          "start_nm":500,"end_nm":800,"points":11,
          "film_thickness_nm":250,"substrate_n":1.5}
    outputs=[]
    for n in (1.,1.33):
        response=client.post("/api/calculate",json={**base,"pore_n":n,"pore_k":0})
        assert response.status_code==200,response.get_json()
        outputs.append(response.get_json())
    dry,wet=outputs
    assert np.max(np.abs(np.array(wet["effective_n"])-dry["effective_n"]))>.005
    assert np.max(np.abs(np.array(wet["R"])-dry["R"]))>1e-6
    for result in outputs:
        np.testing.assert_allclose(np.array(result["R"])+result["T"]+result["A"],1,atol=1e-9)
        assert abs(result["host_fraction"]-.72)<1e-12

def test_pore_dispersion_model():
    response=app.test_client().post("/api/calculate",json={
        "fraction":.05,"pore_fraction":.2,"pore_model":{"model":"cauchy","A":1.33,"B_um2":.002},
        "start_nm":500,"end_nm":800,"points":9})
    assert response.status_code==200,response.get_json()
    result=response.get_json()
    assert result["pore_n"][0]>result["pore_n"][-1]

def test_invalid_three_phase_composition():
    response=app.test_client().post("/api/calculate",json={"fraction":.7,"pore_fraction":.4})
    assert response.status_code==400
