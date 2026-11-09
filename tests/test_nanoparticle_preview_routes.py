"""Research preview route smoke tests: run with flask installed."""
import numpy as np
from nanoparticle_optics.preview_server import app

def test_angular_route_returns_nonempty_curves():
    client=app.test_client()
    r=client.post("/api/angular",json={"wavelength_nm":550,"radius_nm":30,
         "particle_n":2.0,"particle_k":0.1,"host_n":1.5,"distribution":"single"})
    assert r.status_code==200,r.get_json()
    d=r.get_json()
    assert len(d["angle_deg"])==181
    assert len(d["phase_sr_inv"])==181
    assert np.all(np.isfinite(d["phase_sr_inv"]))
    assert -1<=d["g_angular"]<=1

def test_polydisperse_angular_route():
    client=app.test_client()
    r=client.post("/api/angular",json={"wavelength_nm":550,"radius_nm":30,
         "particle_n":2.0,"particle_k":0.1,"host_n":1.5,
         "distribution":"lognormal","spread":0.2})
    assert r.status_code==200,r.get_json()
    assert len(r.get_json()["phase_sr_inv"])==181

def test_ellipsometry_csv_route():
    client=app.test_client()
    csv="wavelength_nm,angle_deg,psi_deg,delta_deg\n"+"\n".join(
        f"{w},{a},30,175" for a in (65,70) for w in (500,550,600,650))
    r=client.post("/api/ellipsometry/parse",json={"csv_text":csv})
    assert r.status_code==200,r.get_json()
    assert len(r.get_json()["angles_deg"])==2
