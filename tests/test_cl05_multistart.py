"""CL-05 regression checks for deterministic causal multistart fitting."""
import numpy as np
import pytest
from scientific_diagnostics import dielectric, fit_causal
from engine import coherent_film

class SyntheticSpectrum:
    def __init__(self, thickness_nm):
        self.wavelength_nm=np.arange(500.,1101.,10.)
        self.thickness_nm=thickness_nm
        self.substrate_n=1.46
        self.ambient_n=1.
        self.angle_deg=0.
        self.polarization='s'
        self.sigma_R=np.ones_like(self.wavelength_nm)
        self.sigma_T=np.ones_like(self.wavelength_nm)
        params=np.array([2.5,.45,1.6,.12,.65,2.2,.18])
        nk=np.sqrt(dielectric(self.wavelength_nm,params,'lorentz',2))
        self.R,self.T=coherent_film(self.wavelength_nm,nk,thickness_nm,1.46,1.,0.,'s')

def test_default_single_start_matches_explicit():
    data=[SyntheticSpectrum(500.),SyntheticSpectrum(1000.)]
    a=fit_causal(data,oscillators=2,max_nfev=100)
    b=fit_causal(data,oscillators=2,max_nfev=100,multistart=1)
    np.testing.assert_allclose(a['result'].x,b['result'].x,rtol=0,atol=1e-12)
    assert len(a['stages']['start_results'])==1

def test_multistart_reproducible_and_not_worse():
    data=[SyntheticSpectrum(500.),SyntheticSpectrum(1000.)]
    single=fit_causal(data,oscillators=2,max_nfev=100,multistart=1)
    a=fit_causal(data,oscillators=2,max_nfev=100,multistart=3,random_seed=42)
    b=fit_causal(data,oscillators=2,max_nfev=100,multistart=3,random_seed=42)
    assert len(a['stages']['start_results'])==3
    assert a['weighted_rmse'] <= single['weighted_rmse']+1e-12
    np.testing.assert_allclose(a['result'].x,b['result'].x,rtol=0,atol=1e-12)

@pytest.mark.parametrize('starts',[0,51])
def test_reject_invalid_starts(starts):
    with pytest.raises(ValueError,match='multistart'):
        fit_causal([SyntheticSpectrum(500.)],multistart=starts)
