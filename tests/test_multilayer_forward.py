"""Multilayer kernel regression tests; run with pytest tests/test_multilayer_forward.py."""
import numpy as np
from multilayer_forward import coherent_stack
from optical_inversion_package.inversion import coherent_film

def test_single_film_matches_existing_forward_model():
    w=np.linspace(400.,1000.,101)
    n=2.1+0.08j+0*w
    for pol in ('s','p'):
        for angle in (0.,35.):
            a,b=coherent_stack(w,[(n,170.)],1.5,1.,angle,pol)
            x,y=coherent_film(w,n,170.,1.5,1.,angle,pol)
            np.testing.assert_allclose(a,x,rtol=1e-11,atol=1e-11)
            np.testing.assert_allclose(b,y,rtol=1e-11,atol=1e-11)

def test_lossless_stack_conserves_power():
    w=np.linspace(420.,980.,113)
    for pol in ('s','p'):
        r,t=coherent_stack(w,[(2.1,125.),(1.8,215.)],1.5,1.,28.,pol)
        np.testing.assert_allclose(r+t,1.,rtol=1e-11,atol=1e-11)

def test_absorbing_stack_passive():
    w=np.linspace(400.,1000.,101)
    r,t=coherent_stack(w,[(2.2+.2j,160.),(1.7,100.)],1.5)
    assert np.all(r>=0) and np.all(t>=0)
    assert np.all(r+t<=1.+1e-10)
