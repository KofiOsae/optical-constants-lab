import pytest
from nanoparticle_optics.meep_colab import MeepJob,make_notebook
from nanoparticle_optics.preview_server import app

@pytest.mark.parametrize("dimension,geometry",[(2,"cylinder"),(3,"sphere")])
def test_notebook_embeds_job(dimension,geometry):
    job=MeepJob(dimension=dimension,geometry=geometry,radius_nm=100)
    notebook=make_notebook(job)
    assert notebook["nbformat"]==4
    src="".join("".join(cell["source"]) for cell in notebook["cells"])
    assert "100.0" in src and "resolution_pixels_per_um" in src
    assert "files.download" in src
    assert "sim.run" in src

def test_invalid_geometry_rejected():
    with pytest.raises(ValueError,match="infinite cylinders"):
        MeepJob(dimension=2,geometry="sphere").validate()

def test_colab_job_requires_user_confirmation():
    client=app.test_client()
    assert client.post("/api/meep/notebook",json={"dimension":2,"geometry":"cylinder"}).status_code==400
    r=client.post("/api/meep/notebook",json={"dimension":3,"geometry":"sphere","radius_nm":100,"user_selected":True})
    assert r.status_code==200,r.get_json()
    assert r.get_json()["notebook"]["nbformat"]==4

def test_unreasonably_large_job_rejected():
    with pytest.raises(ValueError,match="safety cap"):
        MeepJob(dimension=3,geometry="sphere",radius_nm=100,
                padding_um=3,pml_um=.5,resolution_pixels_per_um=200).validate()
