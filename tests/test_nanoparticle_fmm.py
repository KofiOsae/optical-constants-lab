from nanoparticle_optics.materials import MaterialRecord
from nanoparticle_optics.fmm_interchange import export_fmm_material,import_fmm_material
import pytest

def test_roundtrip_and_metadata():
    original=MaterialRecord("test","Synthetic material",[500,600],[1.5,1.6],[.01,.02],
                            composition={"A":1},provenance="synthetic",data_license="original")
    csv,metadata=export_fmm_material(original)
    recovered=import_fmm_material(csv,metadata)
    assert recovered==original

def test_rejects_changed_spectrum():
    original=MaterialRecord("test","Synthetic material",[500,600],[1.5,1.6],[.01,.02],
                            provenance="synthetic",data_license="original")
    csv,metadata=export_fmm_material(original)
    with pytest.raises(ValueError):
        import_fmm_material(csv.replace("1.5","1.9"),metadata)
