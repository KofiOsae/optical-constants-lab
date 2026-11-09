from nanoparticle_optics.preview_server import app
from nanoparticle_optics.assisted_assessment import assess_calculation

def test_assisted_assessment_api():
    r=app.test_client().post("/api/calculate",json={
        "start_nm":500,"end_nm":700,"points":5,"radius_nm":30,
        "host_n":1.33,"host_k":0,"particle_n":1.46,"particle_k":0,
        "fraction":.05,"pore_fraction":.1,"film_thickness_nm":200})
    assert r.status_code==200,r.get_json()
    data=r.get_json()
    a=data["assisted_assessment"]
    assert a["overall"] in ("pass","review")
    names=[c["check"] for c in a["checks"]]
    for name in ("Complex dielectric consistency","Passive effective material","Coherent film energy balance","Mie extinction balance","EMA applicability"):
        assert name in names
    assert all(c["next_step"] for c in a["checks"])

def test_assessment_detects_nonphysical_film():
    d={"wavelength_nm":[500,600],"effective_n":[1.5,1.5],"effective_k":[0,0],
       "effective_epsilon1":[2.25,2.25],"effective_epsilon2":[0,0],
       "R":[.5,.5],"T":[.7,.7],"A":[0,0]}
    result=assess_calculation(d)
    assert result["overall"]=="fail"
