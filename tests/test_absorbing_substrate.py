import numpy as np
from multilayer_forward import coherent_stack

def test_absorbing_substrate_backside_flux():
    w=np.linspace(500,900,31)
    r,t=coherent_stack(w,[(1.8+.06j,100)],2.0+.2j)
    assert np.all(np.isfinite(r)) and np.all(np.isfinite(t))
    assert np.all(r>=0) and np.all(t>=0) and np.all(r+t<=1+1e-9)

def test_finite_lossless_slab_energy():
    w=np.linspace(500,900,31)
    r,t=coherent_stack(w,[],1.5,substrate_thickness_nm=1000,backing_n=1)
    np.testing.assert_allclose(r+t,1,atol=1e-10)

def test_thick_absorbing_slab_has_attenuated_backside():
    w=np.linspace(500,900,31)
    _,t1=coherent_stack(w,[],2+.2j,substrate_thickness_nm=100,backing_n=1)
    _,t2=coherent_stack(w,[],2+.2j,substrate_thickness_nm=1000,backing_n=1)
    assert np.all(t2<t1)

def test_legacy_smooth_substrate_remains_identical():
    w=np.linspace(500,900,31)
    r,t=coherent_stack(w,[(1.8,200)],1.5)
    r2,t2=coherent_stack(w,[(1.8,200)],1.5,substrate_thickness_nm=None)
    np.testing.assert_array_equal(r,r2)
    np.testing.assert_array_equal(t,t2)
