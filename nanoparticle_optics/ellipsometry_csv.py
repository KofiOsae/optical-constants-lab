"""Strict long-form ellipsometry CSV importer for research fitting."""
import csv,io
import numpy as np

COLUMNS=("wavelength_nm","angle_deg","psi_deg","delta_deg")

def parse_ellipsometry_csv(text):
    reader=csv.DictReader(io.StringIO(text.lstrip("\ufeff")))
    if reader.fieldnames is None or any(k not in reader.fieldnames for k in COLUMNS):
        raise ValueError("Required CSV columns: wavelength_nm,angle_deg,psi_deg,delta_deg")
    rows=list(reader)
    if len(rows)<4 or len(rows)>20000:raise ValueError("Expected 4–20000 measurements")
    parsed=[]
    for row in rows:
        try: values=tuple(float(row[k]) for k in COLUMNS)
        except (ValueError,TypeError,KeyError) as e:raise ValueError("Non-numeric or missing ellipsometry measurement") from e
        if not all(np.isfinite(values)):raise ValueError("Non-finite ellipsometry measurement")
        parsed.append(values)
    wl=sorted(set(row[0] for row in parsed))
    angles=sorted(set(row[1] for row in parsed))
    if len(wl)<4 or min(wl)<=0 or min(angles)<=0 or max(angles)>=89:
        raise ValueError("Need at least four positive wavelengths and angles between 0 and 89 degrees")
    psi=np.full((len(angles),len(wl)),np.nan)
    delta=np.full_like(psi,np.nan)
    wi={v:i for i,v in enumerate(wl)};ai={v:i for i,v in enumerate(angles)}
    for w,a,p,d in parsed:
        i=ai[a];j=wi[w]
        if np.isfinite(psi[i,j]):raise ValueError("Duplicate angle/wavelength measurement")
        if not 0<=p<=90:raise ValueError("Psi must be 0–90 degrees")
        psi[i,j]=p;delta[i,j]=d
    if not np.all(np.isfinite(psi)) or not np.all(np.isfinite(delta)):
        raise ValueError("Incomplete angle/wavelength grid; supply a measurement for each combination")
    return {"wavelength_nm":wl,"angles_deg":angles,"psi_deg":psi.tolist(),"delta_deg":delta.tolist()}
