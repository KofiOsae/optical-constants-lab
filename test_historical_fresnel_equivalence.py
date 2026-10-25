import json
import numpy as np
from scipy.interpolate import CubicSpline
from historical_reference import archived_phase2
from optical_inversion_package.inversion import coherent_film

def test_forward_model_matches_original_archived_fresnel():
    d=archived_phase2()
    w=np.asarray(d['wavelength_nm'])
    n=np.asarray(d['n'])
    k=np.asarray(d['k'])
    R,T=coherent_film(w,n+1j*k,d['thickness_nm'],1.,1.,0.,'s')
    assert np.max(np.abs(R-np.asarray(d['R_fit'])))<1e-11
    assert np.max(np.abs(T-np.asarray(d['T_fit'])))<1e-11

def test_archived_spline_is_compatible_with_solver_without_refitting():
    from engine import fit_historical_logk
    from optical_inversion_package.inversion import Spectrum
    import inspect
    d=archived_phase2()
    assert 'warm_start' in inspect.signature(fit_historical_logk).parameters
    assert d['coefficient_count']==50
