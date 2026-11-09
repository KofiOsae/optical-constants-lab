"""Passive transparent Sellmeier and Cauchy material models.

Wavelength is in nanometres externally and micrometres in Sellmeier.
Only transparent n, k=0 is represented; absorption needs a separate model.
"""
import numpy as np

SILICA_MALITSON={"model":"sellmeier","B":[0.6961663,0.4079426,0.8974794],
 "C_um2":[0.00467914826,0.01351206307,97.9340025],
 "min_nm":210,"max_nm":3710,"reference":"Malitson (1965), fused silica; transparent approximation"}

def model_epsilon(model,wavelength_nm):
    wl=np.asarray(wavelength_nm,dtype=float)
    if wl.ndim!=1 or len(wl)<2 or not np.all(np.isfinite(wl)) or np.any(wl<=0):
        raise ValueError("Wavelength array must contain at least two positive finite values")
    kind=model.get("model")
    if kind=="silica_malitson":
        model=SILICA_MALITSON
        kind="sellmeier"
    low=float(model.get("min_nm",0))
    high=float(model.get("max_nm",float("inf")))
    if np.any(wl<low) or np.any(wl>high):
        raise ValueError("Wavelength range exceeds dispersion model validity")
    um=wl/1000.
    if kind=="sellmeier":
        B=np.asarray(model["B"],float);C=np.asarray(model["C_um2"],float)
        if B.ndim!=1 or B.size<1 or B.size>8 or C.shape!=B.shape or not np.all(np.isfinite(B+C)):
            raise ValueError("Sellmeier B and C_um2 arrays must match (1–8 finite terms)")
        denom=um[:,None]**2-C[None,:]
        if np.any(np.abs(denom)<1e-12):raise ValueError("Sellmeier pole in selected range")
        n2=1+np.sum(B[None,:]*um[:,None]**2/denom,axis=1)
        if np.any(n2<=0):raise ValueError("Nonpositive Sellmeier n²; check validity range")
        return n2.astype(complex)
    if kind=="cauchy":
        A=float(model["A"]);B=float(model.get("B_um2",0));C=float(model.get("C_um4",0))
        n=A+B/um**2+C/um**4
        if not np.all(np.isfinite(n)) or np.any(n<=0):raise ValueError("Invalid Cauchy index")
        return (n*n).astype(complex)
    raise ValueError("Unsupported dispersion model")
