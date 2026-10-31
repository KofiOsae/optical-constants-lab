"""End-to-end causal inverse tests using physical R/T forward synthesis."""
import numpy as np
from optical_inversion_package.inversion import Spectrum,coherent_film
from organic_dispersion import organic_nk
from scientific_diagnostics import fit_causal,dielectric,jacobian_report

TRUTH=np.array([2.35,18.,1.5,1.2,.32])
def make_spectrum(thickness,lo=380,hi=950):
    w=np.arange(lo,hi+1,20,dtype=float)
    n,k=organic_nk(w,TRUTH,'tauc_lorentz',1)
    R,T=coherent_film(w,n+1j*k,thickness,1.46)
    return Spectrum(w,R,T,thickness,substrate_n=1.46,sigma_R=.005,sigma_T=.005)

def test_tl_accuracy_near_ground_truth_and_independent_samples():
    specimens=[make_spectrum(180),make_spectrum(300)]
    out=fit_causal(specimens,'tauc_lorentz',1,max_nfev=16,
        parameter_controls={'initial':TRUTH.tolist(),'lower':[2.,12.,1.3,1.,.2],'upper':[2.8,25.,1.7,1.5,.5]})
    assert out['weighted_rmse']<1e-6
    assert abs(out['stages']['coefficients'][2]-1.5)<.01
    w=np.linspace(450,900,50)
    n,k=out['optical_constants'](w)
    rn,rk=organic_nk(w,TRUTH,'tauc_lorentz',1)
    assert np.sqrt(np.mean((n-rn)**2))<.001
    assert np.sqrt(np.mean((k-rk)**2))<.001

def test_tl_reduced_range_does_not_certify_gap():
    specimen=make_spectrum(180,450,690)
    out=fit_causal([specimen],'tauc_lorentz',1,max_nfev=6,
       parameter_controls={'initial':TRUTH.tolist(),'lower':[2.,12.,1.3,1.,.2],'upper':[2.8,25.,1.7,1.5,.5]})
    J=jacobian_report(out['result'])
    assert J['available']
    assert out['weighted_rmse']<1e-6
    # Even an exact fit on a truncated range is not evidence of unique Eg.
    assert max(specimen.wavelength_nm)<1239.841984/TRUTH[2]
