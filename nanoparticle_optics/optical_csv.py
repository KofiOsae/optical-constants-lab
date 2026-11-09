"""Import wavelength-dependent complex optical properties with explicit conventions."""
import csv,io
import numpy as np

def parse_optical_csv(text):
    reader=csv.DictReader(io.StringIO(text.lstrip("\ufeff")))
    names=set(reader.fieldnames or [])
    wave="wavelength_nm"
    mode="nk" if {"n","k"}<=names else "epsilon" if {"epsilon1","epsilon2"}<=names else None
    if wave not in names or mode is None:
        raise ValueError("CSV needs wavelength_nm,n,k or wavelength_nm,epsilon1,epsilon2")
    rows=list(reader)
    if not 2<=len(rows)<=10000:raise ValueError("Require 2–10000 wavelength rows")
    values=[]
    for row in rows:
        try:
            wl=float(row[wave])
            a=float(row["n" if mode=="nk" else "epsilon1"])
            b=float(row["k" if mode=="nk" else "epsilon2"])
        except (ValueError,TypeError,KeyError) as exc:
            raise ValueError("Non-numeric optical properties") from exc
        if not np.all(np.isfinite([wl,a,b])) or wl<=0 or b<0 or (mode=="nk" and a<=0):
            raise ValueError("Invalid wavelength or nonpassive optical constants")
        values.append((wl,a,b))
    values.sort()
    wl=np.array([v[0] for v in values])
    if np.any(np.diff(wl)<=0):raise ValueError("Wavelengths must be unique")
    a=np.array([v[1] for v in values]);b=np.array([v[2] for v in values])
    eps=(a+1j*b)**2 if mode=="nk" else a+1j*b
    return {"wavelength_nm":wl.tolist(),"epsilon1":eps.real.tolist(),"epsilon2":eps.imag.tolist(),"input_mode":mode}

def interpolate_epsilon(record,wavelength_nm):
    x=np.asarray(record["wavelength_nm"],float)
    target=np.asarray(wavelength_nm,float)
    if target.min()<x.min()-1e-8 or target.max()>x.max()+1e-8:
        raise ValueError("Requested wavelengths exceed material CSV coverage; extrapolation disabled")
    return np.interp(target,x,record["epsilon1"])+1j*np.interp(target,x,record["epsilon2"])
