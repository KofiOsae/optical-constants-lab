import numpy as np
from stack_models import ellipsometry
from nanoparticle_optics.ellipsometry_inverse import fit_ellipsometry

def test_synthetic_single_film_recovery():
    wl=np.linspace(450,750,15)
    angles=[60.,70.]
    outputs=[ellipsometry(wl,[{"n":1.8,"k":0.04,"thickness_nm":110.}],
                         substrate_n=1.5,angle_deg=a) for a in angles]
    psi=np.array([o["psi_deg"] for o in outputs])
    delta=np.array([o["delta_deg"] for o in outputs])
    result=fit_ellipsometry(wl,psi,delta,angles,initial=(1.75,.03,105.))
    assert result["success"]
    assert abs(result["n"]-1.8)<1e-3
    assert abs(result["k"]-.04)<1e-3
    assert abs(result["thickness_nm"]-110)<.1
