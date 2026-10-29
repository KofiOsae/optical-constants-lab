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
    print(f'Q01 independent Airy: n_RMSE={n_rmse:.7g}, k_RMSE={k_rmse:.7g}, heldout_RT_RMSE={holdout_rmse:.7g}, branches={fitted["ambiguity"]["distinct_branches"]}')
    assert fitted['result'].success, 'Solver did not converge'
    assert n_rmse<.06, f'Independent n error too large: {n_rmse}'
    assert k_rmse<.025, f'Independent k error too large: {k_rmse}'
    assert holdout_rmse<.012, f'Unseen thickness prediction failed: {holdout_rmse}'
    assert fitted['ambiguity']['near_optimal_branches']>=1
    assert fitted['ambiguity']['status'] in ('competing_solutions','not_established')


def test_independent_lorentz_recovery_withheld_thickness():
    """Blind Lorentz fit against independently specified analytic dielectric data."""
    w=np.linspace(480.,940.,65)
    energy=1239.841984/w
    eps=2.5+.3*(2.1**2)/(2.1**2-energy**2-1j*.2*energy)
    truth=np.sqrt(eps+0j)
    training=[]
    for d in (105.,165.):
        R,T=_independent_airy_rt(w,truth.real,truth.imag,d,substrate_n=1.5)
        training.append((StackSpectrum(w,R,T,(),(),1.5),d))
    fitted=fit_dispersive_stacks(training,model='lorentz',
                                 oscillators=1,max_nfev=250)
    recovered_n,recovered_k=fitted['optical_constants'](w)
    n_error=float(np.sqrt(np.mean((recovered_n-truth.real)**2)))
    k_error=float(np.sqrt(np.mean((recovered_k-truth.imag)**2)))
    R_unseen,T_unseen=_independent_airy_rt(w,truth.real,truth.imag,137.,1.5)
    R_pred,T_pred=coherent_stack(w,[(recovered_n+1j*recovered_k,137.)],1.5)
    holdout_error=float(np.sqrt(np.mean(np.r_[(R_pred-R_unseen),(T_pred-T_unseen)]**2)))
    assert fitted['result'].success
    print(f'Q02 independent Lorentz: n_RMSE={n_error:.7g}, k_RMSE={k_error:.7g}, heldout_RT_RMSE={holdout_error:.7g}')
    assert n_error<.06, f'Lorentz n RMSE {n_error}'
    assert k_error<.025, f'Lorentz k RMSE {k_error}'
    assert holdout_error<.012, f'Lorentz holdout RMSE {holdout_error}'


def test_thin_film_indistinguishable_branches_are_not_reported_unique():
    """Negative control: distinct optical constants have experimentally indistinguishable R/T."""
    w=np.linspace(450.,850.,41)
    aR,aT=_independent_airy_rt(w,1.7,.02,1.,1.5)
    bR,bT=_independent_airy_rt(w,2.6,.01,1.,1.5)
    diff=np.sqrt(np.mean(np.r_[(aR-bR),(aT-bT)]**2))
    assert diff<.001, 'Negative control was not spectrally ambiguous'
    sample=StackSpectrum(w,aR,aT,(),(),1.5,
                         sigma_R=.005,sigma_T=.005)
    fit=fit_dispersive_stacks([(sample,1.)],model='spline',
                              n_basis=4,smoothness=.01,max_nfev=80,
                              optical_bounds=(.1,5.,0.,2.))
    ambiguity=fit['ambiguity']
    # This test is intentionally an identifiability gate rather than a
    # goodness-of-fit gate; a numerically tiny RMSE must not imply unique n,k.
    print(f'Q09 thin-film ambiguity: alt_RT_RMSE={diff:.7g}, distinct_branches={ambiguity["distinct_branches"]}, relative5pct_branches={ambiguity["near_optimal_branches"]}, noise_equivalent_branches={ambiguity["noise_equivalent_branches"]}')
    assert ambiguity['distinct_branches']>=2, ambiguity
    assert ambiguity['noise_equivalent_branches']>=2, ambiguity
    assert ambiguity['status']=='competing_solutions'

def test_previous_fit_warm_start_and_reject_invalid_parameters():
    import pytest
    w=np.linspace(490.,790.,31)
    n=1.92+0.025*(w-640.)/300.
    k=.035+0*w
    r,t=coherent_stack(w,[(n+1j*k,140.)],1.5)
    sample=StackSpectrum(w,r,t,(),(),1.5)
    base=fit_dispersive_stacks([(sample,140.)],model='spline',
                               n_basis=4,max_nfev=75)
    reuse=fit_dispersive_stacks([(sample,140.)],model='spline',
                                n_basis=4,max_nfev=75,
                                initial_parameters=base['dispersion_parameters'])
    assert reuse['warm_start_used'] is True
    assert reuse['multistart']['attempts']>=4
    assert reuse['weighted_rmse']<=base['weighted_rmse']+1.e-5
    assert not base['warm_start_used']
    with pytest.raises(ValueError,match='matching length'):
        fit_dispersive_stacks([(sample,140.)],model='spline',
                              n_basis=4,max_nfev=10,initial_parameters=[1.,2.])
    with pytest.raises(ValueError,match='outside parameter bounds'):
        fit_dispersive_stacks([(sample,140.)],model='spline',
                              n_basis=4,max_nfev=10,
                              initial_parameters=[100.]*8)

def test_progressive_fit_expand_high_wavelength_to_full_spectrum():
    """Stage-one long-band optical fit transfers to a larger B-spline domain."""
    w=np.linspace(500.,900.,41)
    n=1.85+.04*(w-700.)/400.
    k=.025+.005*np.exp(-((w-600.)/85.)**2)
    r,t=coherent_stack(w,[(n+1j*k,180.)],1.5)
    selection=w>=700
    stage_sample=StackSpectrum(w[selection],r[selection],t[selection],(),(),1.5)
    stage1=fit_dispersive_stacks([(stage_sample,180.)],model='spline',
                                n_basis=4,max_nfev=70,optical_bounds=(1.,3.,0.,1.))
    wn=stage_sample.wavelength_nm
    n1,k1=stage1['optical_constants'](wn)
    full_sample=StackSpectrum(w,r,t,(),(),1.5)
    stage2=fit_dispersive_stacks([(full_sample,180.)],model='spline',
                                n_basis=4,max_nfev=100,optical_bounds=(1.,3.,0.,1.),
                                initial_reference={'wavelength_nm':wn.tolist(),
                                                   'n':n1.tolist(),'k':k1.tolist()})
    assert stage2['warm_start_used'] is True
    assert stage2['warm_start_transferred'] is True
    assert stage2['multistart']['attempts']>=4
    nf,kf=stage2['optical_constants'](w)
    rp,tp=coherent_stack(w,[(nf+1j*kf,180.)],1.5)
    err=float(np.sqrt(np.mean(np.r_[rp-r,tp-t]**2)))
    print('PROGRESSIVE_SPLINE full_band_RT_RMSE',err,
          'stage1_success',stage1['result'].success,
          'stage2_success',stage2['result'].success)
    assert stage2['result'].success
    assert err<.01


def test_progressive_transfer_rejects_nonoverlap_or_nonphysical_reference():
    import pytest
    w=np.linspace(500.,800.,31)
    n=np.full_like(w,1.9)
    r,t=coherent_stack(w,[(n+.03j,160.)],1.5)
    sample=StackSpectrum(w,r,t,(),(),1.5)
    with pytest.raises(ValueError,match='overlap'):
        fit_dispersive_stacks([(sample,160.)],model='spline',n_basis=4,max_nfev=10,
                              initial_reference={'wavelength_nm':list(np.arange(900.,1000.,10.)),
                                                 'n':[1.9]*10,'k':[.03]*10})
    with pytest.raises(ValueError,match='valid increasing'):
        fit_dispersive_stacks([(sample,160.)],model='spline',n_basis=4,max_nfev=10,
                              initial_reference={'wavelength_nm':list(np.arange(500.,600.,10.)),
                                                 'n':[1.9]*10,'k':[-.1]*10})
