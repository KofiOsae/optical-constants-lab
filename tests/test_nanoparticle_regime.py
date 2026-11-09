import numpy as np
from nanoparticle_optics.preview_server import app
from nanoparticle_optics.regime import scattering_regime

def test_rayleigh_diagnostic():
    d=scattering_regime(np.array([500.,700.]),5,1.33+0j,np.array([1.46+0j,1.46+0j]))
    assert d["rayleigh_heuristic_all"]
    assert d["solver"]=="independent_homogeneous_sphere_Mie"

def test_aggregate_and_shape_do_not_return_sphere_mie():
    for flags in ({"aggregate":True},{"particle_shape":"spheroid"}):
        response=app.test_client().post("/api/calculate",json={"start_nm":500,"end_nm":700,"points":5,
            "radius_nm":30,"host_n":1.33,"particle_n":1.46,**flags})
        assert response.status_code==200,response.get_json()
        data=response.get_json()
        assert "Csca_nm2" not in data
        assert data["regime_diagnostics"]["solver"]=="external_solver_required"
        assert "external" in data["warning"].lower()

def test_sphere_still_calculates_mie():
    response=app.test_client().post("/api/calculate",json={"start_nm":500,"end_nm":700,"points":5,
        "radius_nm":30,"host_n":1.33,"particle_n":1.46,"particle_shape":"sphere"})
    assert response.status_code==200,response.get_json()
    assert len(response.get_json()["Csca_nm2"])==5
