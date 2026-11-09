import numpy as np
from nanoparticle_optics.preview_server import app

def test_angular_phase_and_polarization_with_silica_coefficients():
    client=app.test_client()
    response=client.post("/api/angular",json={
        "wavelength_nm":550,"radius_nm":40,"host_n":1.33,
        "particle_n":1.46,"particle_k":0,
        "particle_model":{"model":"silica_malitson"},"distribution":"single"})
    assert response.status_code==200,response.get_json()
    d=response.get_json()
    for key in ("angle_deg","phase_sr_inv","perpendicular_nm2_sr","parallel_nm2_sr"):
        assert len(d[key])>=91
        assert np.all(np.isfinite(d[key]))
    assert -1<=d["g_angular"]<=1
    theta=np.deg2rad(d["angle_deg"])
    integral=2*np.pi*np.trapz(np.array(d["phase_sr_inv"])*np.sin(theta),theta)
    np.testing.assert_allclose(integral,1,atol=.01)

def test_angular_rejects_absorbing_host():
    r=app.test_client().post("/api/angular",json={"wavelength_nm":550,"radius_nm":40,
        "particle_n":1.46,"host_material":{"wavelength_nm":[500,600],"epsilon1":[2.2,2.2],"epsilon2":[.1,.1]}})
    assert r.status_code==400
