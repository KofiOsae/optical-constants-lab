import numpy as np
from validation_comparisons_v2 import compare_fits,window_stability,thickness_profile,holdout_metrics

def test_well_constrained_replicates():
    w=np.linspace(450,950,41); n=1.7+0*w;k=.04+0*w
    runs=[dict(parameters=[2.35,18,1.5,1.2,.32],n=n,k=k,weighted_rmse=.0001,thickness_nm=180),
          dict(parameters=[2.36,18.1,1.505,1.195,.325],n=n+.001,k=k+.0005,weighted_rmse=.00011,thickness_nm=185)]
    result=window_stability(runs)
    assert result['status']=='evaluated' and result['max_spectral_deviation_from_mean']<.01
    assert thickness_profile(runs)['thicknesses_nm']==[180,185]
    assert compare_fits(runs,reference={'n':n,'k':k})['reference_rmse_n']<1e-12

def test_ambiguous_equally_good_reconstructions():
    w=np.linspace(450,700,26)
    runs=[dict(parameters=[2,1.5],n=w*0+1.5,k=w*0+.01,weighted_rmse=.001),
          dict(parameters=[3,2.3],n=w*0+2.2,k=w*0+.18,weighted_rmse=.001001)]
    result=compare_fits(runs)
    assert result['max_spectral_deviation_from_mean']>.2
    assert result['rmse_max_over_min']<1.01

def test_held_out_residuals():
    r=holdout_metrics([.1,.2],[.8,.7],[.11,.19],[.79,.71])
    assert r['available'] and r['heldout_points']==4
    assert abs(r['rmse']-.01)<1e-12
