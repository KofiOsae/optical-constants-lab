"""Regression for configurable multilayer multistart on synthetic Lorentz R/T."""
import numpy as np
import pytest
from scientific_diagnostics import dielectric
from multilayer_inverse import StackSpectrum, predict_stack
from multilayer_dispersion import fit_dispersive_stacks

@pytest.fixture(scope="module")
def lorentz_sample():
    wavelength=np.linspace(520.,820.,13)
    parameters=np.array([2.5,0.3,1239.841984/670.,0.2])
    z=np.sqrt(dielectric(wavelength,parameters,"lorentz",1)+0j)
    empty=StackSpectrum(wavelength,None,None,substrate_n=1.5)
    R,T=predict_stack(empty,z.real,np.maximum(z.imag,0.),140.)
    return StackSpectrum(wavelength,np.asarray(R),np.asarray(T),
                         substrate_n=1.5,sigma_R=.02,sigma_T=.02)

def test_multistart_count_best_cost_and_repeatability(lorentz_sample):
    outputs={}
    for count in (1,4,8):
        result=fit_dispersive_stacks([(lorentz_sample,140.)],model="lorentz",
            oscillators=1,n_starts=count,random_seed=2026,max_nfev=35)
        info=result["multistart"]
        assert info["attempts"]==count
        assert len(info["costs"])==count
        assert np.isfinite(info["costs"]).all()
        assert info["selected_start"]==int(np.argmin(info["costs"]))
        assert np.isclose(min(info["costs"]),min(float(v) for v in info["costs"]))
        assert np.isfinite(result["weighted_rmse"])
        outputs[count]=result
    # Adding seeds must never worsen the selected best objective for identical
    # first seeds, even when the numerical fit has not fully converged.
    assert min(outputs[8]["multistart"]["costs"]) <= min(outputs[4]["multistart"]["costs"])+1e-10
    assert min(outputs[4]["multistart"]["costs"]) <= min(outputs[1]["multistart"]["costs"])+1e-10
    replay=fit_dispersive_stacks([(lorentz_sample,140.)],model="lorentz",
        oscillators=1,n_starts=4,random_seed=2026,max_nfev=35)
    np.testing.assert_allclose(replay["multistart"]["costs"],outputs[4]["multistart"]["costs"],
                               rtol=1e-9,atol=1e-9)

@pytest.mark.parametrize("count",[0,21])
def test_invalid_start_count_rejected(lorentz_sample,count):
    with pytest.raises(ValueError,match="Multistart count"):
        fit_dispersive_stacks([(lorentz_sample,140.)],model="lorentz",
            oscillators=1,n_starts=count,max_nfev=10)
