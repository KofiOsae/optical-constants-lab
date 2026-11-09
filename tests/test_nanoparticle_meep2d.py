import json
import pytest
from nanoparticle_optics.meep2d import Meep2DJob, LIMITATION

def test_explicit_selection_required():
    with pytest.raises(PermissionError):
        Meep2DJob(550,40,2).export_json()

def test_2d_job_export():
    data=json.loads(Meep2DJob(550,40,2,user_selected=True).export_json())
    assert data["job"]["dimension"]==2
    assert data["geometry"]=="infinite_cylinder"
    assert "NOT finite 3D spheres" in data["limitations"]

def test_rejects_3d():
    with pytest.raises(ValueError):
        Meep2DJob(550,40,2,dimension=3,user_selected=True).export_json()
