"""Research-only FMM optical material interchange, no direct workbench coupling.

The common contract uses vacuum wavelength in nm and passive n+ik.
Metadata and license are preserved in JSON. CSV is intentionally
restricted to wavelength_nm,n,k; the JSON sidecar preserves provenance.
"""
import csv
import io
import json
import numpy as np
from nanoparticle_optics.materials import MaterialRecord

def export_fmm_material(record):
    record.validate()
    out=io.StringIO()
    writer=csv.writer(out,lineterminator="\n")
    writer.writerow(["wavelength_nm","n","k"])
    writer.writerows(zip(record.wavelength_nm,record.n,record.k))
    sidecar=json.dumps({"schema":"ocl-fmm-optical-material-v1",
                        "material_metadata":json.loads(record.to_json()),
                        "units":{"wavelength":"nm","index":"dimensionless"},
                        "convention":"exp(-i omega t), n+ik"},indent=2)
    return out.getvalue(),sidecar

def import_fmm_material(csv_text,sidecar_text):
    meta=json.loads(sidecar_text)
    if meta.get("schema")!="ocl-fmm-optical-material-v1":
        raise ValueError("Unsupported interchange schema")
    reader=csv.DictReader(io.StringIO(csv_text))
    if reader.fieldnames != ["wavelength_nm","n","k"]:
        raise ValueError("Expected wavelength_nm,n,k columns")
    rows=list(reader)
    if not rows:
        raise ValueError("Empty spectrum")
    record=MaterialRecord.from_json(json.dumps(meta["material_metadata"]))
    wl=np.array([float(r["wavelength_nm"]) for r in rows])
    n=np.array([float(r["n"]) for r in rows])
    k=np.array([float(r["k"]) for r in rows])
    if not (np.array_equal(wl,record.wavelength_nm) and
            np.array_equal(n,record.n) and np.array_equal(k,record.k)):
        raise ValueError("CSV and provenance sidecar disagree")
    return record
