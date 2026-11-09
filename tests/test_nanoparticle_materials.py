import numpy as np
import pytest
from nanoparticle_optics.materials import MaterialRecord

def example():
    return MaterialRecord(material_id="synthetic-zno",name="Synthetic ZnO",
        wavelength_nm=[450,550,650],n=[2.1,2.0,1.95],k=[0.1,0.02,0.001],
        composition={"Zn":1,"O":1},dopants={"Al":"1 at.%"},
        particle_radius_nm=25,particle_radius_std_nm=5,
        synthesis={"method":"hypothetical"},measurement={"technique":"synthetic"},
        provenance="Synthetic demonstration, not measured",data_license="Original synthetic data")

def test_roundtrip():
    original=example()
    loaded=MaterialRecord.from_json(original.to_json())
    assert loaded==original
    assert np.allclose(loaded.permittivity(), (np.array(original.n)+1j*np.array(original.k))**2)

def test_requires_provenance():
    item=example()
    item.provenance=""
    with pytest.raises(ValueError):
        item.validate()

def test_rejects_negative_extinction():
    item=example()
    item.k[0]=-0.01
    with pytest.raises(ValueError):
        item.validate()
