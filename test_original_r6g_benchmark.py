"""Regression checks against archived Workbench export and Phase I fit artifact."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from optical_inversion_package.inversion import coherent_film
ROOT=Path(__file__).resolve().parent

def inputs():
    doc=json.loads((ROOT/'data/original_plot_raw_7.json').read_text())
    curves={v['label']:np.array(v['y'],float) for v in doc['data']['series']}
    w=np.array(doc['data']['series'][0]['x'],float)*1000
    material=pd.read_csv(ROOT/'data/R6G_Singh.csv')
    return doc,w,curves,material

def test_original_geometry():
    doc,w,c,m=inputs()
    assert len(w)==201 and np.isclose(w[0],300) and np.isclose(w[-1],1000)
    assert doc['model']['incident_n']==doc['model']['exit_n']==1
    assert len(doc['model']['layers'])==1
    assert doc['model']['layers'][0]['kind']=='uniform'
    assert doc['model']['layers'][0]['thickness_um']==1

def test_energy_balance():
    _,w,c,_=inputs()
    assert np.max(np.abs(c['R']+c['T']+c['A']-1))<1e-12

def test_independent_forward_model():
    _,w,c,m=inputs()
    wm=m['wavelength_um'].to_numpy()*1000
    nf=np.interp(w,wm,m['n'])+1j*np.interp(w,wm,m['k'])
    r,t=coherent_film(w,nf,1000)
    assert np.sqrt(np.mean((r-c['R'])**2))<0.005
    assert np.sqrt(np.mean((t-c['T'])**2))<0.005

def test_phase1_archived_fit():
    _,w,c,_=inputs()
    x=np.load(ROOT/'phase1_fit_verified.npz')
    assert abs(float(x['q'][-1])*1000-999.38)<0.02
    # Stored fit inputs match the archived spectrum, not a newly generated replacement.
    assert np.max(np.abs(np.interp(x['wl']*1000,w,c['R'])-x['R']))<1e-9
    assert np.max(np.abs(np.interp(x['wl']*1000,w,c['T'])-x['T']))<1e-9
