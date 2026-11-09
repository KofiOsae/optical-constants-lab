import numpy as np
import pytest
from nanoparticle_optics.distributions import radius_bins

@pytest.mark.parametrize("kind",["single","normal","lognormal","gamma","weibull","bimodal"])
def test_distribution_weights(kind):
    r,w=radius_bins(kind,30,.2,count=31)
    assert np.all(r>0) and np.all(w>=0)
    assert abs(sum(w)-1)<1e-10

def test_histogram_number_weights():
    r,w=radius_bins("histogram",30,histogram_radii_nm=[10,20],histogram_number_weights=[1,3])
    assert w.tolist()==[.25,.75]
