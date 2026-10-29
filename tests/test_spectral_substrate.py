import numpy as np
import pytest
from multilayer_api import substrate_index
from multilayer_forward import coherent_stack

def test_spectral_complex_substrate_interpolated():
    wavelength=np.array([500.,600.,700.])
    config={'substrate_n':1.5,'substrate_k':0.,
            'substrate_ri':{'wavelength_nm':[400.,800.],
                            'n':[1.7,2.1],'k':[.01,.09]}}
    ns=substrate_index(config,wavelength)
    np.testing.assert_allclose(ns.real,[1.8,1.9,2.0])
    np.testing.assert_allclose(ns.imag,[.03,.05,.07])
    R,T=coherent_stack(wavelength,[(1.8+.03j,130.)],substrate_n=ns)
    assert np.all(np.isfinite(R)) and np.all(np.isfinite(T))
    assert np.all(R>=0) and np.all(T>=0)
    assert np.all(R+T<=1+1e-10)

def test_substrate_RI_requires_coverage():
    with pytest.raises(ValueError,match='cover'):
        substrate_index({'substrate_ri':{'wavelength_nm':[550.,650.],
             'n':[1.5,1.5],'k':[0.,0.]}},np.array([500.,600.,700.]))
