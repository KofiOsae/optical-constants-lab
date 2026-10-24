from benchmark_api import run_original_inverse

def test_reference_assisted_roundtrip():
    r=run_original_inverse(950)
    assert abs(r['reference_assisted_thickness_nm']-1000)<.01
    assert r['reference_assisted_spectral_rmse']<1e-10

def test_blind_fit_initialization_sensitivity():
    near=run_original_inverse(1000)
    far=run_original_inverse(950)
    assert near['blind_spectral_rmse']<1e-4
    assert far['blind_spectral_rmse']>near['blind_spectral_rmse']
