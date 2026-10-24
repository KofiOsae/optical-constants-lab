import numpy as np
from engine import FitSettings, make_problem, harris_hawks
from optical_inversion_package.inversion import Spectrum, coherent_film

def test_swanepoel_seed_is_accepted_and_clipped():
    w=np.linspace(450,850,60)
    n=np.full_like(w,2.0)
    R,T=coherent_film(w,n+0j,800,1.5)
    sample=Spectrum(w,R,T,800,substrate_n=1.5,name='synthetic')
    settings=FitSettings(n_basis=6,fit_thickness=True,initial_n=(w,n),initial_thickness_nm=805)
    residual,x0,lo,hi,_,_=make_problem([sample],settings)
    assert np.all(x0>=lo) and np.all(x0<=hi)
    assert np.isclose(x0[-1],805)
    assert np.isfinite(residual(x0)).all()

def test_hho_bounded_and_reproducible():
    low=np.array([-2.,-2.]); high=np.array([2.,2.])
    objective=lambda x: float(np.sum((x-.25)**2))
    a,ha,ea=harris_hawks(objective,low,high,population=8,iterations=5,seed=7)
    b,hb,eb=harris_hawks(objective,low,high,population=8,iterations=5,seed=7)
    assert np.allclose(a,b)
    assert np.allclose(ha,hb)
    assert ea==eb
    assert np.all(a>=low) and np.all(a<=high)
    assert np.all(np.diff(ha)<=1e-12)

def test_swanepoel_seeded_fresnel_spline_recovers_synthetic_film():
    from engine import fit
    w=np.linspace(420,880,81)
    truth_n=1.9; truth_k=0.015; d=780.
    R,T=coherent_film(w,np.full_like(w,truth_n+1j*truth_k,dtype=complex),d,1.5)
    sample=Spectrum(w,R,T,d,substrate_n=1.5,name='known-film')
    settings=FitSettings(n_basis=6,fit_thickness=True,initial_n=(w,np.full_like(w,truth_n)),
                         initial_thickness_nm=d,optimizer='least_squares',max_nfev=120)
    result=fit([sample],settings)
    n,k=result['optical_constants'](w)
    assert result['result'].success
    assert abs(result['thicknesses'][0]-d)<2.
    assert np.mean(abs(n-truth_n))<.01
    assert np.mean(abs(k-truth_k))<.01
    assert result['weighted_rmse']<.05
