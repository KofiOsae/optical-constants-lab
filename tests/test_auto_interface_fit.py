import numpy as np
from optical_inversion_package.inversion import Spectrum,coherent_film
from scientific_diagnostics import fit_causal,dielectric
from multilayer_inverse import StackSpectrum,predict_stack
from multilayer_dispersion import fit_dispersive_stacks

def _test_data():
    w=np.linspace(450,850,51)
    coefficients=[2.5,.8,2.3,.12]
    z=np.sqrt(dielectric(w,coefficients,'lorentz',1)+0j)
    return w,z.real,z.imag,coefficients

def _bounded_coefficients(c):
    c=np.asarray(c,float)
    delta=np.array([.002,.002,.002,.002])
    return {'initial':c.tolist(),'lower':(c-delta).tolist(),'upper':(c+delta).tolist()}

def test_causal_shared_interface_mixing_recovery():
    w,n,k,p=_test_data()
    samples=[]
    for thickness in (140.,280.):
        R,T=coherent_film(w,n+1j*k,thickness,1.5,roughness_nm=2.)
        samples.append(Spectrum(w,R,T,thickness,substrate_n=1.5))
    out=fit_causal(samples,model='lorentz',oscillators=1,max_nfev=90,
                   parameter_controls=_bounded_coefficients(p),
                   fit_interface_mixing=True,interface_mixing_bounds=(0.,5.))
    assert out['result'].success
    assert np.isfinite(out['weighted_rmse'])
    assert out['weighted_rmse']<1.e-3
    assert abs(out['fitted_interface_mixing_nm']-2.)<.5
    assert len(out['stages']['coefficients'])==4

def test_multilayer_shared_interface_mixing_recovery():
    w,n,k,p=_test_data()
    samples=[]
    for thickness in (140.,280.):
        generated=StackSpectrum(w,None,None,substrate_n=1.5,roughness_nm=(2.,2.))
        R,T=predict_stack(generated,n,k,thickness)
        sample=StackSpectrum(w,R,T,substrate_n=1.5)
        samples.append((sample,thickness))
    out=fit_dispersive_stacks(samples,model='lorentz',oscillators=1,max_nfev=90,
                              parameter_controls=_bounded_coefficients(p),
                              fit_interface_mixing=True,interface_mixing_bounds=(0.,5.))
    assert out['result'].success
    assert out['weighted_rmse']<1.e-3
    assert abs(out['fitted_interface_mixing_nm']-2.)<.5
    assert len(out['dispersion_parameters'])==4
