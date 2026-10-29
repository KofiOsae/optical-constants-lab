"""Scientific checks for the shared organic dielectric foundations."""
import numpy as np
import pytest
from organic_dispersion import organic_dielectric,organic_nk,MODEL_PARAMETER_COUNTS

@pytest.mark.parametrize("model,p",[
 ("modified_lorentz",[2.5,0.8,2.3,0.15,0.4]),
 ("kk_gaussian",[2.5,0.8,2.3,0.15]),
 ("voigt",[2.5,0.8,2.3,0.12,0.15]),
 ("vibronic",[2.5,0.8,2.1,0.18,0.12,1.0]),
])
def test_organic_passivity_and_finite_nk(model,p):
    wavelength=np.linspace(350.,1100.,151)
    eps=organic_dielectric(wavelength,p,model)
    n,k=organic_nk(wavelength,p,model)
    assert np.isfinite(eps).all()
    assert np.all(eps.imag>=-1e-12)
    assert np.all(n>0)
    assert np.all(k>=-1e-12)
    np.testing.assert_allclose((n+1j*k)**2,eps,rtol=1e-10,atol=1e-10)

def test_zero_strength_is_constant_background():
    for model,size in MODEL_PARAMETER_COUNTS.items():
        params={"modified_lorentz":[0,2,.2,.3],"kk_gaussian":[0,2,.2],
                "voigt":[0,2,.1,.2],"vibronic":[0,2,.2,.1,.5]}[model]
        assert len(params)==size
        eps=organic_dielectric(np.array([400.,800.]),[2.25,*params],model)
        np.testing.assert_allclose(eps,2.25)

def test_invalid_negative_strength_rejected():
    with pytest.raises(ValueError):
        organic_dielectric([500.],[2.5,-1,2.3,.2,.4],"modified_lorentz")
