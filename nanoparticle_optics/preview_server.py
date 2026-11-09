"""Standalone research preview: python -m nanoparticle_optics.preview_server

This server is NOT registered in production server.py and has no database
or external job execution. For local research use only.
"""
from flask import Flask, request, jsonify, send_from_directory
from pathlib import Path
import numpy as np
from nanoparticle_optics.effective_medium import maxwell_garnett, bruggeman
from nanoparticle_optics.mie import mie_sphere
from nanoparticle_optics.angular import mie_anisotropy, ensemble_anisotropy
from nanoparticle_optics.phase_function import sphere_phase_function
from nanoparticle_optics.ensemble import ensemble_cross_sections
from nanoparticle_optics.distributions import radius_bins
from nanoparticle_optics.materials import MaterialRecord
from nanoparticle_optics.fmm_interchange import export_fmm_material
from nanoparticle_optics.meep2d import Meep2DJob
from nanoparticle_optics.ellipsometry_inverse import fit_ellipsometry
from nanoparticle_optics.ellipsometry_csv import parse_ellipsometry_csv

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
        result={"wavelength_nm":wl.tolist(),"volume_fraction":f,"effective_n":[float(nk.real)]*len(wl),
                "effective_k":[float(nk.imag)]*len(wl),"effective_epsilon1":[float(eps.real)]*len(wl),
                "effective_epsilon2":[float(eps.imag)]*len(wl),"method":method}
        # A transparent nonabsorbing host is required for the Mie implementation.
        if host.imag==0:
            radius=float(p.get("radius_nm",30))
            if p.get("distribution","single")!="single":
                radii,weights=radius_bins(p.get("distribution","lognormal"),radius,float(p.get("spread",.2)),count=61,second_radius_nm=p.get("second_radius_nm"),mixture_fraction=float(p.get("mixture_fraction",.5)),histogram_radii_nm=p.get("histogram_radii_nm"),histogram_number_weights=p.get("histogram_number_weights"))
                m=ensemble_cross_sections(wl,radii,weights,inc,host.real)
                result.update({k:m[k].tolist() for k in ("Cext_nm2","Csca_nm2","Cabs_nm2")})
                result["radii_nm"]=radii.tolist()
                result["number_weights"]=weights.tolist()
                result["volume_weights"]=m["volume_weights"].tolist()
                result["mean_geometric_area_nm2"]=float(np.pi*np.dot(weights,radii**2))
                angular=ensemble_anisotropy(wl,radii,weights,inc,host.real,number_density_nm3=f/m["mean_particle_volume_nm3"])
                result.update({k:angular[k].tolist() for k in ("g","mu_s_nm_inv","mu_s_prime_nm_inv")})
            else:
                m=mie_sphere(wl,radius,inc,host.real)
                result.update({k:m[k].tolist() for k in ("Cext_nm2","Csca_nm2","Cabs_nm2")})
                angular=mie_anisotropy(wl,radius,inc,host.real)
                density=f/(4*np.pi/3*radius**3)
                result["g"]=angular["g"].tolist()
                result["mu_s_nm_inv"]=(density*m["Csca_nm2"]).tolist()
                result["mu_s_prime_nm_inv"]=(density*m["Csca_nm2"]*(1-angular["g"])).tolist()
                result["mean_geometric_area_nm2"]=float(np.pi*radius**2)
        else:
            result["warning"]="Mie scattering disabled: host must be lossless."
        return jsonify(result)
    except (ValueError,TypeError,OverflowError,ZeroDivisionError) as exc:
        return jsonify({"error":str(exc)}),400

@app.post("/api/material/validate")
def validate_material():
    try:
        record=MaterialRecord.from_json(__import__("json").dumps(request.get_json(force=True)))
        csv,sidecar=export_fmm_material(record)
        return jsonify({"valid":True,"material_id":record.material_id,"fmm_csv":csv,"fmm_metadata_json":sidecar})
    except (ValueError,TypeError,KeyError) as exc:
        return jsonify({"error":str(exc)}),400

@app.post("/api/meep2d/export")
def export_meep2d():
    try:
        job=Meep2DJob(**request.get_json(force=True))
        return jsonify(__import__("json").loads(job.export_json()))
    except (ValueError,TypeError,PermissionError) as exc:
        return jsonify({"error":str(exc)}),400

@app.post("/api/ellipsometry/parse")
def ellipsometry_parse_api():
    try:
        p=request.get_json(force=True)
        if len(p["csv_text"])>2000000:raise ValueError("CSV exceeds 2 MB")
        return jsonify(parse_ellipsometry_csv(p["csv_text"]))
    except (ValueError,TypeError,KeyError) as exc:
        return jsonify({"error":str(exc)}),400

@app.post("/api/ellipsometry/fit")
def fit_ellipsometry_api():
    try:
        p=request.get_json(force=True)
        return jsonify(fit_ellipsometry(p["wavelength_nm"],p["psi_deg"],p["delta_deg"],p["angles_deg"],
                 substrate_n=float(p.get("substrate_n",1.5)),
                 initial=p.get("initial",[1.7,.02,100]),
                 lower=p.get("lower",[1.01,0,1]),upper=p.get("upper",[5,3,3000])))
    except (ValueError,TypeError,KeyError) as exc:
        return jsonify({"error":str(exc)}),400

@app.post("/api/distribution")
def distribution_api():
    try:
        p=request.get_json(force=True)
        radii,weights=radius_bins(p.get("distribution","single"),float(p.get("radius_nm",30)),
          float(p.get("spread",.2)),count=61,second_radius_nm=p.get("second_radius_nm"),
          mixture_fraction=float(p.get("mixture_fraction",.5)),
          histogram_radii_nm=p.get("histogram_radii_nm"),histogram_number_weights=p.get("histogram_number_weights"))
        volume=weights*radii**3
        return jsonify({"radii_nm":radii.tolist(),"number_weights":weights.tolist(),
          "volume_weights":(volume/volume.sum()).tolist()})
    except (ValueError,TypeError,KeyError,OverflowError) as exc:
        return jsonify({"error":str(exc)}),400

@app.post("/api/angular")
def angular_api():
    try:
        p=request.get_json(force=True)
        wl=float(p["wavelength_nm"])
        host=float(p.get("host_n",1.5))
        inc=complex(float(p["particle_n"]),float(p.get("particle_k",0)))
        radius=float(p["radius_nm"])
        if p.get("distribution","single")=="single":
            out=sphere_phase_function(wl,radius,inc,host)
            return jsonify({k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in out.items()})
        radii,weights=radius_bins(p["distribution"],radius,float(p.get("spread",.2)),count=61,
              second_radius_nm=p.get("second_radius_nm"),mixture_fraction=float(p.get("mixture_fraction",.5)),
              histogram_radii_nm=p.get("histogram_radii_nm"),histogram_number_weights=p.get("histogram_number_weights"))
        mean=None
        for r,w in zip(radii,weights):
            if w==0:continue
            v=sphere_phase_function(wl,r,inc,host)
            if mean is None:
                mean={k:np.zeros_like(v[k]) for k in ("perpendicular_nm2_sr","parallel_nm2_sr")}
            for key in mean:mean[key]+=w*v[key]
        theta=v["angle_deg"]
        total=.5*(mean["perpendicular_nm2_sr"]+mean["parallel_nm2_sr"])
        csca=2*np.pi*np.trapz(total*np.sin(np.deg2rad(theta)),np.deg2rad(theta))
        phase=total/csca
        return jsonify({"angle_deg":theta.tolist(),"phase_sr_inv":phase.tolist(),
          "perpendicular_nm2_sr":mean["perpendicular_nm2_sr"].tolist(),
          "parallel_nm2_sr":mean["parallel_nm2_sr"].tolist(),
          "Csca_angular_nm2":float(csca),
          "g_angular":float(2*np.pi*np.trapz(phase*np.cos(np.deg2rad(theta))*np.sin(np.deg2rad(theta)),np.deg2rad(theta)))})
    except (ValueError,TypeError,KeyError,OverflowError,ZeroDivisionError) as exc:
        return jsonify({"error":str(exc)}),400

if __name__=="__main__":
    app.run(host="127.0.0.1",port=8766,debug=False)
