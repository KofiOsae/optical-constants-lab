import numpy as np
import pytest
from multilayer_forward import coherent_stack
from optical_inversion_package.inversion import coherent_film

def test_smooth_film_exact_legacy():
    w=np.linspace(450,950,41)
    film=1.75+.12j
    base=coherent_film(w,film,150.,1.5)
    zero=coherent_film(w,film,150.,1.5,roughness_nm=0.)
    mixed=coherent_stack(w,[(film,150.)],1.5,roughness_nm=[0.,0.])
    for a,b,c in zip(base,zero,mixed):
        np.testing.assert_allclose(a,b,rtol=0,atol=0)
        np.testing.assert_allclose(a,c,rtol=1e-12,atol=1e-12)

def test_mixing_layer_has_effect_and_is_passive():
    w=np.linspace(450,950,51)
    base=coherent_stack(w,[(1.8+.08j,140.)],1.5)
    mixed=coherent_stack(w,[(1.8+.08j,140.)],1.5,roughness_nm=[4.,5.])
    assert max(np.max(np.abs(a-b)) for a,b in zip(base,mixed))>1e-6
    r,t=mixed
    assert np.all(r>=0) and np.all(t>=0)
    assert np.all(r+t<=1+1e-10)

def test_interface_count_and_negative_thickness_rejected():
    w=np.linspace(450,950,11)
    with pytest.raises(ValueError):
        coherent_stack(w,[(1.8,140.)],1.5,roughness_nm=[2.])
    with pytest.raises(ValueError):
        coherent_stack(w,[(1.8,140.)],1.5,roughness_nm=[0.,-1.])
