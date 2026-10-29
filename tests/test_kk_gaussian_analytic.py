import numpy as np
from scipy.integrate import quad
from organic_dispersion import _analytic_gaussian_lorentz, organic_dielectric

def test_broad_gaussian_against_adaptive_integration():
    E=np.array([1.6,2.0,2.3,2.45,2.8])
    strength=.15
    center=2.3
    gamma=.005
    sigma=.15
    analytic=_analytic_gaussian_lorentz(E,strength,center,gamma,sigma)
    for energy,actual in zip(E,analytic):
        def f(c):
            weight=np.exp(-.5*((c-center)/sigma)**2)/(sigma*np.sqrt(2*np.pi))
            return weight*strength*c*c/(c*c-energy**2-1j*gamma*energy)
        expected=complex(quad(lambda c:f(c).real,-np.inf,np.inf,epsabs=1e-9)[0],
                         quad(lambda c:f(c).imag,-np.inf,np.inf,epsabs=1e-9)[0])
        np.testing.assert_allclose(actual,expected,atol=3e-7,rtol=2e-7)
        assert actual.imag>=0

def test_narrow_gaussian_continuity_and_passivity():
    w=np.linspace(400,1000,301)
    params=[2.5,.15,2.3,.002]
    eps=organic_dielectric(w,params,'kk_gaussian')
    assert np.all(np.isfinite(eps))
    assert np.all(eps.imag>=0)
    assert np.max(eps.imag)>0
