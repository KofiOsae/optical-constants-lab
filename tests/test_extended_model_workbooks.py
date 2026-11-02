"""Numerical model workbook must reconstruct fitted constants from parameters."""
import numpy as np
import pytest
from openpyxl import load_workbook
from model_workbook import build_numerical_model_workbook
from organic_dispersion import organic_nk

@pytest.mark.parametrize('model,params',[
 ('modified_lorentz',[2.5,.4,2.1,.15,.3]),
 ('kk_gaussian',[2.5,.4,2.1,.15]),
 ('voigt',[2.5,.4,2.1,.12,.15]),
 ('vibronic',[2.5,.4,2.1,.18,.12,1.]),
 ('tauc_lorentz',[2.5,.4,1.4,.5,.2]),
])
def test_organic_export_matches_exact_model(model,params):
    wavelength=np.array([550.,650.,850.])
    book=load_workbook(build_numerical_model_workbook(model,params,wavelength))
    ws=book['Reconstruction']
    n,k=organic_nk(wavelength,params,model,1)
    np.testing.assert_allclose([ws.cell(row,4).value for row in range(2,5)],n,rtol=1e-10)
    np.testing.assert_allclose([ws.cell(row,5).value for row in range(2,5)],k,rtol=1e-10)
    assert 'not live Excel' in book['External_reconstruction']['A2'].value

def test_spline_export_domain_and_coefficients():
    w=np.array([500.,650.,800.])
    out=build_numerical_model_workbook('spline',[1.7,1.8,1.9,2.,.01,.02,.03,.04],w,{'domain_min_nm':500,'domain_max_nm':800})
    book=load_workbook(out)
    assert book['Reconstruction']['D2'].value==pytest.approx(1.7)
    assert book['Reconstruction']['D4'].value==pytest.approx(2.)
    assert book['Parameters']['C9'].value==pytest.approx(.04)
