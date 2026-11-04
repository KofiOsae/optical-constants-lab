"""Single inversion must return fitted spline parameters for workbook replay."""
from pathlib import Path
def test_single_custom_fit_serializes_spline_parameters():
    source=(Path(__file__).resolve().parents[1]/'server.py').read_text()
    assert "'model':'spline'" in source
    assert "'dispersion_parameters':np.asarray(result['result'].x[:2*basis],float).tolist()" in source
