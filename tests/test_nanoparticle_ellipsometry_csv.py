import pytest
from nanoparticle_optics.ellipsometry_csv import parse_ellipsometry_csv

def test_multiangle_import():
    rows=["wavelength_nm,angle_deg,psi_deg,delta_deg"]
    rows += [f"{w},{a},30,175" for a in [65,70] for w in [500,550,600,650]]
    data=parse_ellipsometry_csv("\n".join(rows))
    assert data["angles_deg"]==[65.,70.]
    assert len(data["psi_deg"])==2
    assert len(data["psi_deg"][0])==4

def test_missing_grid_rejected():
    rows=["wavelength_nm,angle_deg,psi_deg,delta_deg"]
    rows += [f"{w},{a},30,175" for a in [65,70] for w in [500,550,600,650] if not(a==70 and w==600)]
    with pytest.raises(ValueError,match="Incomplete"):parse_ellipsometry_csv("\n".join(rows))

def test_duplicate_rejected():
    rows=["wavelength_nm,angle_deg,psi_deg,delta_deg"]
    rows += [f"{w},65,30,175" for w in [500,550,600,650]]
    rows += ["500,65,30,175"]
    with pytest.raises(ValueError,match="Duplicate"):parse_ellipsometry_csv("\n".join(rows))
