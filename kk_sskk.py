"""Finite-band KK and singly-subtractive KK diagnostics.

Energy E=hc/lambda in eV. The output is a band-limited principal-value
estimate, NOT a full-spectrum causal reconstruction. Unknown tails remain
uncontrolled. SSKK requires an independently known real-epsilon anchor.
"""
import numpy as np

HC_EV_NM=1239.841984

def _grid(wavelength_nm, epsilon2):
    w=np.asarray(wavelength_nm,float)
    e2=np.asarray(epsilon2,float)
    if w.ndim!=1 or w.size<12 or e2.shape!=w.shape or not np.all(np.isfinite(w)) or not np.all(np.isfinite(e2)) or np.any(w<=0):
        raise ValueError('Require >=12 finite positive wavelengths and aligned epsilon2')
    energy=HC_EV_NM/w
    order=np.argsort(energy)
    energy=energy[order];e2=e2[order]
    if np.any(np.diff(energy)<=0):raise ValueError('Wavelengths must be unique')
    return energy,e2

def _pv_integral(energy,epsilon2,queries):
    """PV integral of E*epsilon2(E)/(E^2-query^2) over finite band."""
    a,b=energy[0],energy[-1]
    g=energy*epsilon2
    values=[]
    for q in np.asarray(queries,float):
        if not a<q<b:raise ValueError('Query/anchor must lie strictly inside measured energy band')
        gq=np.interp(q,energy,g)
        denominator=energy*energy-q*q
        numerator=g-gq
        # Analytic limit at the removable singularity (interpolated slope).
        slope=np.gradient(g,energy)
        nearest=np.argmin(np.abs(energy-q))
        safe=np.abs(denominator)>1e-12*max(1,q*q)
        regular=np.empty_like(energy)
        regular[safe]=numerator[safe]/denominator[safe]
        regular[~safe]=np.interp(q,energy,slope)/(2*q)
        analytic=gq/(2*q)*np.log(((b-q)*(a+q))/((b+q)*(q-a)))
        values.append(np.trapezoid(regular,energy)+analytic)
    return np.asarray(values)

def kk_finite_band(wavelength_nm,epsilon2,query_wavelength_nm=None,epsilon_inf=1.):
    """Finite-band epsilon1 estimate; assumes zero epsilon2 outside band."""
    energy,e2=_grid(wavelength_nm,epsilon2)
    if query_wavelength_nm is None:
        queries=energy[1:-1]
    else:
        queries=HC_EV_NM/np.atleast_1d(np.asarray(query_wavelength_nm,float))
    values=float(epsilon_inf)+(2/np.pi)*_pv_integral(energy,e2,queries)
    return {'wavelength_nm':(HC_EV_NM/queries).tolist(),
            'epsilon1_estimate':values.tolist(),
            'warning':'Finite-band KK depends on unknown spectral tails and epsilon_inf; not a causality proof.'}

def sskk_finite_band(wavelength_nm,epsilon2,anchor_wavelength_nm,anchor_epsilon1,
                     query_wavelength_nm=None):
    """Finite-band singly subtractive estimate with independent epsilon1 anchor."""
    energy,e2=_grid(wavelength_nm,epsilon2)
    anchor=HC_EV_NM/float(anchor_wavelength_nm)
    queries=energy[1:-1] if query_wavelength_nm is None else HC_EV_NM/np.atleast_1d(np.asarray(query_wavelength_nm,float))
    integrals=_pv_integral(energy,e2,np.r_[queries,anchor])
    values=float(anchor_epsilon1)+(2/np.pi)*(integrals[:-1]-integrals[-1])
    return {'wavelength_nm':(HC_EV_NM/queries).tolist(),
            'epsilon1_estimate':values.tolist(),
            'anchor_wavelength_nm':float(anchor_wavelength_nm),
            'anchor_epsilon1':float(anchor_epsilon1),
            'warning':'SSKK reduces some finite-band sensitivity but does not eliminate unknown tails; anchor must be independent.'}
