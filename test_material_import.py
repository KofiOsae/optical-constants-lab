import pytest
from material_import import inspect_material

def test_nk_upload_without_reference():
    result = inspect_material('unknown.csv', b'wavelength_nm,n,k\n500,1.5,0.01\n600,1.6,0.02')
    assert result['kind'] == 'optical_constants'
    assert result['columns']['k'] == [0.01, 0.02]

def test_rt_um():
    result = inspect_material('unknown.csv', b'wavelength_um,R,T\n0.5,0.2,0.7\n0.6,0.3,0.6')
    assert result['kind'] == 'film_spectrum'
    assert result['wavelength_nm'] == [500., 600.]

def test_ellipsometry_preview():
    result = inspect_material('unknown.csv', b'wavelength_nm,psi,delta\n500,45,120\n600,46,121')
    assert result['kind'] == 'ellipsometry'
    assert result['warnings']

def test_ambiguous_units_rejected():
    with pytest.raises(ValueError, match='Ambiguous wavelength'):
        inspect_material('unknown.csv', b'wavelength,n,k\n500,1.5,0\n600,1.6,0')
