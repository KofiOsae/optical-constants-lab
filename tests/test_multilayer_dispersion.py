import numpy as np
from multilayer_inverse import StackSpectrum
from multilayer_dispersion import fit_dispersive_stacks
from multilayer_forward import coherent_stack
from scientific_diagnostics import dielectric

def test_joint_spline_two_stacks():
    w=np.linspace(470.,870.,65)
    n=2.05+0.06*(w-650.)/400.
    k=.055+0*w
    samples=[]
    for d in (110.,180.):
        before=((1.7,85.),)
        after=((1.9,55.),)
        r,t=coherent_stack(w,[*before,(n+1j*k,d),*after],1.52)
        samples.append((StackSpectrum(w,r,t,before,after,1.52),d))
    out=fit_dispersive_stacks(samples,model='spline',n_basis=6,max_nfev=180)
    nn,kk=out['optical_constants'](w)
    assert out['result'].success
    assert np.sqrt(np.mean((nn-n)**2))<.03
    assert np.sqrt(np.mean((kk-k)**2))<.02
    assert 'rank' in out['sensitivity']

def test_joint_lorentz_known_answer():
    w=np.linspace(480.,940.,80)
    coeff=[2.5,.3,2.1,.2]
    z=np.sqrt(dielectric(w,coeff,'lorentz',1)+0j)
    samples=[]
    for d in (105.,165.):
        before=((1.8,75.),)
        r,t=coherent_stack(w,[*before,(z,d)],1.5)
        samples.append((StackSpectrum(w,r,t,before,(),1.5),d))
    out=fit_dispersive_stacks(samples,model='lorentz',oscillators=1,max_nfev=250)
    assert out['result'].success
    assert out['weighted_rmse']<.03


def test_multilayer_branch_diagnostic_is_reference_independent():
    """Branch diagnostics must use fitted spectra only, with no reference n,k."""
    w=np.linspace(500.,800.,21)
    n=np.full_like(w,1.8)
    k=np.full_like(w,.08)
    r,t=coherent_stack(w,[(n+1j*k,120.)],1.5)
    sample=StackSpectrum(w,r,t,(),(),1.5)
    out=fit_dispersive_stacks([(sample,120.)],model='spline',
                              n_basis=4,max_nfev=70,
                              optical_bounds=(.01,8.,0.,10.))
    report=out['ambiguity']
    assert report['distinct_branches']>=1
    assert report['near_optimal_branches']>=1
    assert report['status'] in ('competing_solutions','not_established')
    assert len(report['probe_wavelength_nm'])==51
    assert all(len(b['n_probe'])==51 and len(b['k_probe'])==51
               for b in report['branches'])
    assert all(b['relative_cost_increase']>=-1e-8 for b in report['branches'])
    assert 'reference' not in out


def _independent_airy_rt(w, n, k, thickness_nm, substrate_n=1.5):
    """Independent normal-incidence Airy amplitude formula; not coherent_stack."""
    N=np.asarray(n)+1j*np.asarray(k)
    r01=(1-N)/(1+N)
    r12=(N-substrate_n)/(N+substrate_n)
    t01=2/(1+N)
    t12=2*N/(N+substrate_n)
    phase=np.exp(2j*np.pi*N*thickness_nm/np.asarray(w))
    denominator=1+r01*r12*phase**2
    R=np.abs((r01+r12*phase**2)/denominator)**2
    T=substrate_n*np.abs(t01*t12*phase/denominator)**2
    return R,T


def test_independent_airy_blind_recovery_and_unseen_thickness():
    """No use of truth/reference RI during fitting or branch selection."""
    w=np.linspace(470.,850.,65)
    truth_n=1.85+0.045*(w-650.)/380.
    truth_k=.055+.012*np.exp(-((w-590.)/95.)**2)
    rng=np.random.default_rng(20261029)
    training=[]
    for d in (115.,205.):
        R,T=_independent_airy_rt(w,truth_n,truth_k,d)
        R=np.clip(R+rng.normal(0,.0005,len(w)),0,1)
        T=np.clip(T+rng.normal(0,.0005,len(w)),0,1)
        training.append((StackSpectrum(w,R,T,(),(),1.5,sigma_R=.0005,sigma_T=.0005),d))
    fitted=fit_dispersive_stacks(training,model='spline',n_basis=5,
                                 smoothness=.01,max_nfev=140,
                                 optical_bounds=(1.,3.,0.,.5))
    recovered_n,recovered_k=fitted['optical_constants'](w)
    n_rmse=float(np.sqrt(np.mean((recovered_n-truth_n)**2)))
    k_rmse=float(np.sqrt(np.mean((recovered_k-truth_k)**2)))
    # Unseen specimen is generated independently, and is never passed to fit.
    truth_R,truth_T=_independent_airy_rt(w,truth_n,truth_k,157.)
    predicted_R,predicted_T=coherent_stack(w,[(recovered_n+1j*recovered_k,157.)],1.5)
    holdout_rmse=float(np.sqrt(np.mean(np.r_[(predicted_R-truth_R),
                                               (predicted_T-truth_T)]**2)))
    assert fitted['result'].success, 'Solver did not converge'
    assert n_rmse<.06, f'Independent n error too large: {n_rmse}'
    assert k_rmse<.025, f'Independent k error too large: {k_rmse}'
    assert holdout_rmse<.012, f'Unseen thickness prediction failed: {holdout_rmse}'
    assert fitted['ambiguity']['near_optimal_branches']>=1
    assert fitted['ambiguity']['status'] in ('competing_solutions','not_established')
