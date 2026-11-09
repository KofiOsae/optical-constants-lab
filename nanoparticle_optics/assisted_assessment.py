"""Assisted engineering assessment; heuristic flags are not literature validation."""
import numpy as np

def assess_calculation(data):
    checks=[]
    def add(name,status,detail,action):
        checks.append({"check":name,"status":status,"detail":detail,"next_step":action})
    wl=np.asarray(data["wavelength_nm"],float)
    nk=np.asarray(data["effective_n"],float)+1j*np.asarray(data["effective_k"],float)
    eps=np.asarray(data["effective_epsilon1"],float)+1j*np.asarray(data["effective_epsilon2"],float)
    residual=float(np.max(np.abs(nk*nk-eps)))
    add("Complex dielectric consistency","pass" if residual<1e-7 else "fail",
        f"max |(n+ik)^2-epsilon|={residual:.3g}","Check square-root branch and units.")
    passive=bool(np.all(np.isfinite(eps)) and np.min(eps.imag)>=-1e-8)
    add("Passive effective material","pass" if passive else "fail",
        f"minimum Im(epsilon)={np.min(eps.imag):.3g}","Inspect constituent k, oscillator damping, mixing model and branch.")
    if all(k in data for k in ("R","T","A")):
        r,t,a=(np.asarray(data[k],float) for k in ("R","T","A"))
        energy=float(np.max(np.abs(r+t+a-1)))
        physical=bool(np.all(r>=-1e-7) and np.all(t>=-1e-7) and np.all(a>=-1e-7))
        add("Coherent film energy balance","pass" if energy<1e-6 and physical else "fail",
            f"max |R+T+A-1|={energy:.3g}; nonnegative={physical}",
            "Check TMM normalization, substrate and absorption convention. Does not assess diffuse scattering.")
    if all(k in data for k in ("Cext_nm2","Csca_nm2","Cabs_nm2")):
        ext,sca,absb=(np.asarray(data[k],float) for k in ("Cext_nm2","Csca_nm2","Cabs_nm2"))
        err=float(np.max(np.abs(ext-sca-absb))/max(1,float(np.max(np.abs(ext)))))
        ok=bool(err<1e-6 and np.min(absb)>=-1e-6 and np.min(sca)>=-1e-6)
        add("Mie extinction balance","pass" if ok else "fail",
            f"relative |Cext-Csca-Cabs|={err:.3g}","Check sphere solver, sign convention and numerical convergence.")
    if "g" in data:
        g=np.asarray(data["g"],float)
        add("Scattering anisotropy","pass" if np.all(np.isfinite(g)) and np.max(np.abs(g))<=1+1e-6 else "fail",
            f"g range {np.min(g):.3g} to {np.max(g):.3g}",
            "Check angular integration and scattering weighting across distributions.")
    frac=float(data.get("volume_fraction",0))
    pore=float(data.get("pore_fraction",0))
    add("EMA applicability","review" if frac+pore>.15 else "pass",
        f"particle fraction={frac:.3g}, pore fraction={pore:.3g}",
        "At higher loading or connected/aggregated phases, compare alternative mixing models and microstructure measurements.")
    if data.get("film_model")=="coherent_specular_TMM_no_diffuse_scattering":
        add("Specular vs diffuse losses","review","Film R/T/A excludes diffuse scattering and haze.",
            "Measure integrating-sphere total/specular R/T and assess whether independent-sphere scattering is negligible.")
    if data.get("regime_diagnostics"):
        for warning in data["regime_diagnostics"]["warnings"]:
            add("Particle regime","review",warning,"Check size distribution, geometry and solver assumptions.")
    if any(c["status"]=="fail" for c in checks):overall="fail"
    elif any(c["status"]=="review" for c in checks):overall="review"
    else:overall="pass"
    return {"overall":overall,"checks":checks,"disclaimer":"Engineering self-checks only; pass does not establish experimental accuracy or published validation."}
