import numpy as np
from organic_dispersion import organic_dielectric, organic_parameter_setup

def test_tauc_lorentz_passivity_and_subgap():
    wavelength=np.linspace(300,1100,151)
    p=np.array([2.4,15.0,1.6,1.1,0.35])
    eps=organic_dielectric(wavelength,p,'tauc_lorentz',1)
    E=1239.841984/wavelength
    assert np.all(np.isfinite(eps))
    assert np.all(eps.imag>=-1e-12)
    assert np.max(np.abs(eps.imag[E<=1.6]))<1e-12
    assert np.max(eps.imag[E>1.6])>0
    x,lo,hi=organic_parameter_setup('tauc_lorentz',1,300,1100)
    assert len(x)==5 and np.all(lo<=x) and np.all(x<=hi)

def test_tauc_lorentz_two_oscillator_additivity():
    w=np.linspace(400,900,55)
    b=2.2
    q1=[8,1.2,1.4,.25]
    q2=[4,1.8,.8,.3]
    together=organic_dielectric(w,[b,*q1,*q2],'tauc_lorentz',2)
    single1=organic_dielectric(w,[b,*q1],'tauc_lorentz',1)
    single2=organic_dielectric(w,[b,*q2],'tauc_lorentz',1)
    assert np.allclose(together,single1+single2-b,rtol=1e-9,atol=1e-9)
