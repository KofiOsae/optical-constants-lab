import numpy as np
import pytest
from material_import import inspect_material
from custom_spectra import parse_spectrum

def test_preview_and_fit_use_same_units_and_sorted_data():
    raw = b'wavelength_um,R_percent,T_percent\n' + b''.join(f'{x:.3f},20,70\n'.encode() for x in np.linspace(.4,.8,9))
    preview=inspect_material('sample.csv',raw)
    spectrum=parse_spectrum('sample.csv',raw,500,1.5)
    np.testing.assert_allclose(preview['wavelength_nm'],spectrum.wavelength_nm)
    np.testing.assert_allclose(preview['columns']['R'],spectrum.R)
    np.testing.assert_allclose(spectrum.T,.7)

def test_ambiguous_units_fail_consistently():
    raw=b'wavelength,R\n'+b''.join(f'{x},0.2\n'.encode() for x in range(400,490,10))
    with pytest.raises(ValueError,match='Ambiguous'):
        inspect_material('sample.csv',raw)
    with pytest.raises(ValueError,match='Ambiguous'):
        parse_spectrum('sample.csv',raw,500)
    assert len(parse_spectrum('sample.csv',raw,500,wavelength_unit='nm').wavelength_nm)==9

def test_nk_preview_does_not_require_r6g_or_fit():
    raw=b'wavelength_nm,n,k\n'+b''.join(f'{x},1.5,0.01\n'.encode() for x in range(400,490,10))
    preview=inspect_material('material.csv',raw)
    assert preview['kind']=='optical_constants'
    with pytest.raises(ValueError,match='requires reflectance'):
        parse_spectrum('material.csv',raw,500)
