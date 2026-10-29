import numpy as np
from kk_sskk import kk_finite_band, sskk_finite_band

def test_finite_band_kk_and_sskk_numpy_2():
    w=np.linspace(400.,1100.,201)
    e=1239.841984/w
    eps2=.15*np.exp(-((e-2.3)/.2)**2)
    kk=kk_finite_band(w,eps2,epsilon_inf=2.)
    assert len(kk['wavelength_nm'])==199
    assert np.all(np.isfinite(kk['epsilon1_estimate']))
    sskk=sskk_finite_band(w,eps2,anchor_wavelength_nm=700.,anchor_epsilon1=2.8)
    assert len(sskk['wavelength_nm'])==199
    assert np.all(np.isfinite(sskk['epsilon1_estimate']))
