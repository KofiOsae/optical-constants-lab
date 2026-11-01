"""Excel reconstruction regression for shared causal models (Single/Joint/Multilayer)."""
from io import BytesIO
from openpyxl import load_workbook
import pytest
from model_workbook import build_model_workbook

@pytest.mark.parametrize('model,parameters',[
    ('lorentz',[2.5,.4,2.1,.15]),
    ('drude_lorentz',[2.5,.3,.15,.4,2.1,.15]),
])
def test_causal_excel_has_editable_parameters_and_formulas(model,parameters):
    file=build_model_workbook(model,parameters,[500.,650.,800.],{'workflow':'multilayer'})
    book=load_workbook(BytesIO(file.getvalue()))
    assert list(book.sheetnames)==['Parameters','Reconstruction','Metadata','Instructions']
    assert book['Parameters']['B2'].value==pytest.approx(2.5)
    ws=book['Reconstruction']
    assert ws['A2'].value==500.
    assert ws['B2'].value.startswith('=')
    for col in 'CDEF':
        assert ws[f'{col}2'].data_type=='f'
    assert book['Metadata']['B2'].value==model

def test_unsupported_model_never_claims_formula_reconstruction():
    with pytest.raises(ValueError,match='requires Lorentz'):
        build_model_workbook('tauc_lorentz',[2.5,.4,1.,.5,.2],[600.])
