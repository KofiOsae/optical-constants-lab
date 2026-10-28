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
