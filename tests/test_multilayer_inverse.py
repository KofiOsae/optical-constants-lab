import numpy as np
from multilayer_forward import coherent_stack
from multilayer_inverse import StackSpectrum,fit_constant_unknown

def test_known_layers_recover_constant_unknown():
    wl=np.linspace(460.,900.,89)
    before=((1.75,90.),)
    after=((2.05,70.),)
    n,k,d=2.35,.11,145.
    r,t=coherent_stack(wl,[*before,(n+1j*k,d),*after],1.52)
    sample=StackSpectrum(wl,r,t,before,after,1.52)
    fit=fit_constant_unknown(sample,d,n0=2.3,k0=.08)
    assert fit['result'].success
    assert abs(fit['n']-n)<1e-5
    assert abs(fit['k']-k)<1e-5
    assert fit['weighted_rmse']<1e-7

def test_unknown_thickness_recovery():
    wl=np.linspace(450.,950.,101)
    before=((1.8,100.),)
    n,k,d=2.25,.09,130.
    r,t=coherent_stack(wl,[*before,(n+1j*k,d)],1.5)
    sample=StackSpectrum(wl,r,t,before,(),1.5)
    fit=fit_constant_unknown(sample,125.,n0=2.2,k0=.07,
                             fit_thickness=True,thickness_bounds=(115.,145.))
    assert fit['result'].success
    assert abs(fit['n']-n)<1e-4
    assert abs(fit['k']-k)<1e-4
    assert abs(fit['thickness_nm']-d)<.01

def test_three_finite_layers_unknown_middle_thickness_and_index():
    """Known A / unknown B / known C with simultaneous n,k,thickness recovery."""
    wl=np.linspace(460.,900.,89)
    known_before=((1.75,90.),)
    known_after=((2.05,70.),)
    n,k,d=2.35,.11,145.
    rr,tt=coherent_stack(wl,[*known_before,(n+1j*k,d),*known_after],1.52)
    spec=StackSpectrum(wl,rr,tt,known_before,known_after,1.52)
    out=fit_constant_unknown(spec,142.,n0=2.3,k0=.08,
                             fit_thickness=True,thickness_bounds=(130.,160.))
    print('THREE_LAYER_INVERSION',out['n'],out['k'],out['thickness_nm'],out['weighted_rmse'])
    assert out['result'].success
    assert abs(out['n']-n)<1e-4
    assert abs(out['k']-k)<1e-4
    assert abs(out['thickness_nm']-d)<.01
