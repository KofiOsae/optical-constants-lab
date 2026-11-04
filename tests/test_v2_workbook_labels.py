from openpyxl import load_workbook
from model_workbook import build_numerical_model_workbook

def test_tauc_lorentz_parameters_have_physical_names_and_units():
    coefficients = [2.35, 18., 1.5, 1.2, 0.32]
    book = load_workbook(build_numerical_model_workbook('tauc_lorentz', coefficients, [450., 550., 690.]))
    params = book['Parameters']
    assert [params.cell(i, 2).value for i in range(2, 7)] == [
        'Background ε∞', 'Oscillator 1 · Amplitude A', 'Oscillator 1 · Optical gap Eg',
        'Oscillator 1 · E0−Eg', 'Oscillator 1 · Broadening C']
    assert [params.cell(i, 3).value for i in range(2, 7)] == coefficients
    assert [params.cell(i, 4).value for i in range(2, 7)] == [
        'dimensionless', 'model amplitude', 'eV', 'eV', 'eV']
    assert book['Reconstruction'].max_row == 4

def test_spline_export_labels_n_and_k_coefficients_separately():
    coefficients = [1.6] * 4 + [0.] * 4
    book = load_workbook(build_numerical_model_workbook(
        'spline', coefficients, [450., 550., 690.],
        {'domain_min_nm': 450., 'domain_max_nm': 690.}))
    params = book['Parameters']
    assert params.cell(2, 2).value == 'n spline basis 1'
    assert params.cell(5, 2).value == 'n spline basis 4'
    assert params.cell(6, 2).value == 'k spline basis 1'
    assert params.cell(9, 2).value == 'k spline basis 4'
