"""Qualified classical analysis; no claim of unique n,k recovery."""
import numpy as np
from scipy.signal import find_peaks
from scipy.interpolate import PchipInterpolator

def swanepoel(wavelength_nm, transmission, substrate_n):
    w=np.asarray(wavelength_nm,float); t=np.asarray(transmission,float); s=float(substrate_n)
    if len(w)<12 or np.any(np.diff(w)<=0) or not (1<s<5): raise ValueError("Swanepoel needs ordered wavelengths and a known transparent substrate n>1")
    if np.any((t<=0)|(t>1)) or np.any(~np.isfinite(t)): raise ValueError("Transmission must be fractions in (0,1]")
    prom=max(.005,.08*float(np.ptp(t))); dist=max(2,len(w)//35)
    maxima,_=find_peaks(t,prominence=prom,distance=dist)
    minima,_=find_peaks(-t,prominence=prom,distance=dist)
    if min(len(maxima),len(minima))<2: raise ValueError("At least two resolved maxima and minima required")
    mask=(w>=max(w[maxima[0]],w[minima[0]]))&(w<=min(w[maxima[-1]],w[minima[-1]]))
    x=w[mask]
    if len(x)<4: raise ValueError("Fringe envelopes have insufficient overlap")
    upper=PchipInterpolator(w[maxima],t[maxima])(x)
    lower=PchipInterpolator(w[minima],t[minima])(x)
    good=(upper>lower+1e-5)&(lower>0)&(upper<=1.01)
    x=x[good];upper=upper[good];lower=lower[good]
    if len(x)<4: raise ValueError("Fringe envelopes cross or contrast is insufficient")
    N=2*s*(upper-lower)/(upper*lower)+(s*s+1)/2
    if np.any(N*N<s*s): raise ValueError("Negative envelope discriminant")
    n=np.sqrt(N+np.sqrt(N*N-s*s))
    return {"method":"Swanepoel envelopes","wavelength_nm":x.tolist(),"n":n.tolist(),"T_max":upper.tolist(),"T_min":lower.tolist(),"note":"Normal-incidence weakly absorbing uniform film on a thick transparent substrate; envelope-derived n only, not k. Check substrate back-face assumptions."}

def tauc(wavelength_nm,absorbance=None,transmission=None,thickness_nm=None,exponent=2,fit_range_ev=None):
    w=np.asarray(wavelength_nm,float)
    d=float(thickness_nm)
    if d<=0: raise ValueError("Positive film thickness required")
    if absorbance is None:
        if transmission is None: raise ValueError("Absorbance or transmission required")
        T=np.asarray(transmission,float)
        if np.any((T<=0)|(T>1)): raise ValueError("Transmission must be fractions in (0,1]")
        A=-np.log10(T); assumption="Absorbance estimated from T only; Fresnel losses and interference neglected."
    else:
        A=np.asarray(absorbance,float);assumption="Measured decadic absorbance; scattering/reflection corrections not verified."
    if len(w)!=len(A) or np.any(w<=0) or np.any(~np.isfinite(A)) or np.any(A<0): raise ValueError("Invalid spectral data")
    E=1239.841984/w; alpha=np.log(10)*A/(d*1e-7)
    y=(alpha*E)**float(exponent)
    order=np.argsort(E);E=E[order];y=y[order]
    out={"energy_ev":E.tolist(),"tauc_y":y.tolist(),"note":assumption+" Optical gap depends on transition model and selected linear interval."}
    if fit_range_ev is not None:
        lo,hi=map(float,fit_range_ev);sel=(E>=lo)&(E<=hi)
        if sel.sum()<4: raise ValueError("Choose a Tauc fit interval with at least four points")
        m,b=np.polyfit(E[sel],y[sel],1)
        if m<=0: raise ValueError("Tauc slope is nonpositive; select the absorption edge")
        fit=m*E[sel]+b
        out.update(bandgap_ev=float(-b/m),r2=float(1-np.sum((y[sel]-fit)**2)/max(np.sum((y[sel]-y[sel].mean())**2),1e-30)),fit_range_ev=[lo,hi])
    return out
