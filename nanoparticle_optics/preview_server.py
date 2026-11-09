"""Standalone research preview: python -m nanoparticle_optics.preview_server

This server is NOT registered in production server.py and has no database
or external job execution. For local research use only.
"""
from flask import Flask, request, jsonify, send_from_directory
from pathlib import Path
import numpy as np
from nanoparticle_optics.effective_medium import maxwell_garnett, bruggeman
from nanoparticle_optics.mie import mie_sphere
from nanoparticle_optics.ensemble import ensemble_cross_sections, lognormal_radius_bins

app=Flask(__name__)
ROOT=Path(__file__).resolve().parent

@app.get("/")
def index():
    return send_from_directory(ROOT,"preview.html")

@app.post("/api/calculate")
def calculate():
    try:
        p=request.get_json(force=True)
        wl=np.linspace(float(p.get("start_nm",400)),float(p.get("end_nm",900)),int(p.get("points",101)))
        if len(wl)>1001 or len(wl)<2 or wl[0]<=0 or wl[-1]<=wl[0]:
            raise ValueError("Invalid wavelength interval")
        host=complex(float(p.get("host_n",1.5)),float(p.get("host_k",0)))
        inc=complex(float(p.get("particle_n",2)),float(p.get("particle_k",0.1)))
        f=float(p.get("fraction",.1))
        method=p.get("method","maxwell_garnett")
        if method not in ("maxwell_garnett","bruggeman"):
            raise ValueError("Unknown effective-medium model")
        eps=(maxwell_garnett if method=="maxwell_garnett" else bruggeman)(host**2,inc**2,f)
        nk=np.sqrt(eps+0j)
        if nk.imag<0: nk=-nk
        result={"wavelength_nm":wl.tolist(),"effective_n":[float(nk.real)]*len(wl),
                "effective_k":[float(nk.imag)]*len(wl),"method":method}
        # A transparent nonabsorbing host is required for the Mie implementation.
        if host.imag==0:
            radius=float(p.get("radius_nm",30))
            if p.get("distribution","single")=="lognormal":
                radii,weights=lognormal_radius_bins(radius,float(p.get("geometric_std",1.2)),count=81)
                m=ensemble_cross_sections(wl,radii,weights,inc,host.real)
                result.update({k:m[k].tolist() for k in ("Cext_nm2","Csca_nm2","Cabs_nm2")})
            else:
                m=mie_sphere(wl,radius,inc,host.real)
                result.update({k:m[k].tolist() for k in ("Cext_nm2","Csca_nm2","Cabs_nm2")})
        else:
            result["warning"]="Mie scattering disabled: host must be lossless."
        return jsonify(result)
    except (ValueError,TypeError,OverflowError,ZeroDivisionError) as exc:
        return jsonify({"error":str(exc)}),400

if __name__=="__main__":
    app.run(host="127.0.0.1",port=8766,debug=False)
