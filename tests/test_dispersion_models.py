import numpy as np
import pytest
from dispersion_models import cauchy,sellmeier,fit_cauchy

def test_cauchy_exact_recovery():
    w=np.linspace(450,1000,111)
    n=cauchy(w,1.45,.006,.0002)
    fit=fit_cauchy(w,n,3)
    assert np.allclose(fit['coefficients'],[1.45,.006,.0002],rtol=0,atol=1e-11)
    assert fit['rmse']<1e-12

def test_cauchy_two_parameter():
    w=np.linspace(500,900,41)
    assert np.allclose(fit_cauchy(w,cauchy(w,1.5,.01),2)['coefficients'],[1.5,.01])

def test_sellmeier_fused_silica():
    # Malitson fused silica coefficients; wavelength in micrometres, C in um^2.
    B=[.6961663,.4079426,.8974794]
    C=[.0684043**2,.1162414**2,9.896161**2]
    p=[v for pair in zip(B,C) for v in pair]
    assert abs(float(sellmeier([587.6],p)[0])-1.45846)<2e-4

def test_sellmeier_rejects_poles():
    with pytest.raises(ValueError):
        sellmeier([500],[1,.25,0,1,0,2])
