import numpy as np
import pytest
from nanoparticle_optics.optical_csv import parse_optical_csv,interpolate_epsilon
from nanoparticle_optics.preview_server import app

def test_nk_and_epsilon_import():
    nk=parse_optical_csv("wavelength_nm,n,k\n500,1.5,0\n700,1.6,0.1")
    eps=parse_optical_csv("wavelength_nm,epsilon1,epsilon2\n500,2.25,0\n700,2.55,0.32")
    np.testing.assert_allclose(interpolate_epsilon(nk,[500]),[2.25+0j])
    np.testing.assert_allclose(interpolate_epsilon(eps,[700]),[2.55+.32j])
    with pytest.raises(ValueError,match="extrapolation"):
        interpolate_epsilon(nk,[400,600])

def test_dispersive_effective_epsilon_and_spectral_mie():
    wl=[500,600,700]
    host=parse_optical_csv("wavelength_nm,n,k\n500,1.4,0\n700,1.5,0")
    particle=parse_optical_csv("wavelength_nm,epsilon1,epsilon2\n500,3.5,0.1\n700,4.5,0.3")
    r=app.test_client().post("/api/calculate",json={"start_nm":500,"end_nm":700,
      "points":3,"fraction":0.05,"host_material":host,"particle_material":particle,
      "radius_nm":30,"distribution":"single"})
    assert r.status_code==200,r.get_json()
    d=r.get_json()
    assert len(d["Csca_nm2"])==3
    assert d["effective_epsilon1"][0]!=d["effective_epsilon1"][-1]
    np.testing.assert_allclose((np.array(d["effective_n"])+1j*np.array(d["effective_k"]))**2,
        np.array(d["effective_epsilon1"])+1j*np.array(d["effective_epsilon2"]),rtol=1e-10)
