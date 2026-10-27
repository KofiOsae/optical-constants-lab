"""Reproducible Lorentz known-answer benchmark. Run: python -m tests.lorentz_benchmark
This is a diagnostic benchmark, not an assertion of successful recovery.
"""
import json
import numpy as np
from optical_inversion_package.inversion import Spectrum, coherent_film
from scientific_diagnostics import dielectric, fit_causal, jacobian_report

def run():
    wavelengths=np.linspace(400.,1000.,151)
    # eps_inf; (strength, resonance energy eV, damping eV) x 2
    truth=np.array([2.4,1.2,2.1,.30,.7,3.0,.48])
    eps=dielectric(wavelengths,truth,'lorentz',2)
    index=np.sqrt(eps+0j)
    spectra=[]
    for thickness in (100.,300.):
        R,T=coherent_film(wavelengths,index,thickness,1.50)
        spectra.append(Spectrum(wavelength_nm=wavelengths,R=R,T=T,
            thickness_nm=thickness,substrate_n=1.50,
            sigma_R=.005,sigma_T=.005,name=f'lorentz_{thickness:.0f}nm'))
    result=fit_causal(spectra,'lorentz',2,max_nfev=400,
                      fit_thickness=False,multistart=5)
    nfit,kfit=result['optical_constants'](wavelengths)
    ntrue,ktrue=index.real,index.imag
    channel_errors={}
    for spec in spectra:
        Rfit,Tfit=coherent_film(wavelengths,nfit+1j*kfit,spec.thickness_nm,1.50)
        channel_errors[spec.name]={
            'R_rmse':float(np.sqrt(np.mean((Rfit-spec.R)**2))),
            'T_rmse':float(np.sqrt(np.mean((Tfit-spec.T)**2)))}
    report={'model':'two-oscillator Lorentz','points':len(wavelengths),
        'thickness_nm':[100,300],'true_coefficients':truth.tolist(),
        'fitted_coefficients':result['stages']['coefficients'],
        'weighted_RT_rmse':result['weighted_rmse'],
        'n_rmse':float(np.sqrt(np.mean((nfit-ntrue)**2))),
        'k_rmse':float(np.sqrt(np.mean((kfit-ktrue)**2))),
        'channels':channel_errors,
        'jacobian':jacobian_report(result['result']),
        'multistart':result['stages']['start_results']}
    print(json.dumps(report,indent=2))
    return report

if __name__=='__main__':
    run()
