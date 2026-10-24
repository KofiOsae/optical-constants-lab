import numpy as np
from stack_models import stack_spectrum,ellipsometry
from optical_inversion_package.inversion import coherent_film

def test_zero_layers_fresnel():
    w=np.linspace(400,800,9)
    r=stack_spectrum(w,[],substrate_n=1.5)
    assert np.allclose(r['R'],.04,atol=1e-10)
    assert np.allclose(r['T'],.96,atol=1e-10)

def test_single_layer_matches_existing_solver():
    w=np.linspace(450,850,20)
    layer=[{'n':1.9,'k':.04,'thickness_nm':120}]
    r=stack_spectrum(w,layer,substrate_n=1.5)
    R,T=coherent_film(w,1.9+.04j,120,substrate_n=1.5)
    assert np.allclose(r['R'],R,atol=1e-9)
    assert np.allclose(r['T'],T,atol=1e-9)

def test_passive_energy_balance():
    w=np.linspace(400,900,31)
    r=stack_spectrum(w,[{'n':2.0,'k':.1,'thickness_nm':200}],substrate_n=1.5)
    assert np.min(r['A'])>=-1e-9
    assert np.allclose(r['R']+r['T']+r['A'],1)

def test_ellipsometry_finite():
    w=np.linspace(400,900,20)
    d=ellipsometry(w,[{'n':1.8,'k':.02,'thickness_nm':100}])
    assert np.all(np.isfinite(d['psi_deg']))
    assert np.all(np.isfinite(d['delta_deg']))
