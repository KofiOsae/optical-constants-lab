import numpy as np
import pytest
from nanoparticle_optics.effective_medium import maxwell_garnett, bruggeman

@pytest.mark.parametrize("model", [maxwell_garnett, bruggeman])
def test_endpoints(model):
    assert model(2+0.1j, 5+0.2j, 0) == 2+0.1j
    assert model(2+0.1j, 5+0.2j, 1) == 5+0.2j

@pytest.mark.parametrize("model", [maxwell_garnett, bruggeman])
def test_identical_media(model):
    assert np.allclose(model(2+0.1j, 2+0.1j, [0,.2,.8,1]), 2+0.1j)

def test_mg_dilute_limit():
    h, i, f = 2.0, 4.0, 1e-6
    expected = h + 3*f*h*(i-h)/(i+2*h)
    assert abs(maxwell_garnett(h,i,f)-expected) < 1e-11

def test_bruggeman_residual():
    h, i, f = 2+0.1j, 5+0.2j, .3
    e = bruggeman(h,i,f)
    residual = f*(i-e)/(i+2*e)+(1-f)*(h-e)/(h+2*e)
    assert abs(residual) < 1e-12

@pytest.mark.parametrize("model", [maxwell_garnett, bruggeman])
def test_invalid_fraction(model):
    with pytest.raises(ValueError):
        model(2,4,1.1)
