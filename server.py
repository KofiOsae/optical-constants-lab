from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from flask import Flask, jsonify, request, send_from_directory
import os
import numpy as np
import pandas as pd
from engine import FitSettings,fit,simulate_four,compare_prediction
from optical_inversion_package.inversion import Spectrum,coherent_film
app=Flask(__name__,static_folder='web',static_url_path='')
app.config['MAX_CONTENT_LENGTH']=2_000_000
DATA=ROOT/'data'/'R6G_Singh.csv'

def reference():
    df=pd.read_csv(DATA)
    cols={c.lower().strip():c for c in df.columns}
    def col(names):
        for name in names:
            if name in cols:return np.asarray(df[cols[name]],dtype=float)
        for name in names:
            if len(name)>2:
                for key,c in cols.items():
                    if name in key:return np.asarray(df[c],dtype=float)
        raise ValueError('Missing column '+str(names))
    w=col(['wavelength_um','wavelength_nm','wavelength','wl','lambda']);n=col(['n']);k=col(['k'])
    if 'wavelength_um' in cols or (('wavelength_nm' not in cols) and np.nanmax(w)<10):w=w*1000
    order=np.argsort(w)
    if len(w)<3 or not np.all(np.isfinite(w)) or not np.all(np.isfinite(n)) or not np.all(np.isfinite(k)) or np.any(np.diff(w[order])<=0):raise ValueError('Invalid R6G reference dataset')
    return w[order],n[order],k[order]

def samples():
    w,n,k=reference()
    mask=(w>=400)&(w<=950)
    w,n,k=w[mask],n[mask],k[mask]
    idx=np.unique(np.linspace(0,len(w)-1,min(61,len(w)),dtype=int))
    return simulate_four(w[idx],n[idx],k[idx]),(w[idx],n[idx],k[idx])

from benchmark_api import run_original_inverse
from custom_spectra import parse_spectrum
from material_import import inspect_material
from analysis_methods import swanepoel, tauc
from stack_models import stack_spectrum, ellipsometry

@app.post('/api/original-inverse')
def original_inverse():
    try:
        payload=request.get_json(silent=True) or {}
        return jsonify(run_original_inverse(float(payload.get('initial_thickness_nm',950))))
    except Exception as exc:
        return jsonify({'error':str(exc)}),400

@app.errorhandler(413)
def upload_too_large(error):
    return jsonify({'error':'Upload too large. Maximum file size is 1 MB.'}),413

@app.get('/health')
def health():return jsonify({'status':'ok','service':'optical-constants-lab','version':'0.3.9'})

@app.get('/')
def home():return send_from_directory(ROOT/'web','index.html')

@app.post('/api/import-preview')
def import_preview():
    try:
        uploaded = request.files.get('spectrum')
        if not uploaded:
            raise ValueError('Choose a CSV, TSV, TXT or JSON file')
        result = inspect_material(
            uploaded.filename, uploaded.read(),
            request.form.get('wavelength_unit', 'auto'),
            request.form.get('signal_unit', 'auto'))
        return jsonify(result)
    except (ValueError, TypeError, KeyError) as exc:
        return jsonify({'error': str(exc)}), 400

@app.post('/api/classical-analysis')
def classical_analysis():
    try:
        uploaded=request.files.get('spectrum')
        if not uploaded: raise ValueError('Select a spectrum file')
        preview=inspect_material(uploaded.filename,uploaded.read(),request.form.get('wavelength_unit','auto'),request.form.get('signal_unit','auto'))
        material={'wavelength_nm':preview['wavelength_nm'],**preview['columns']}
        method=request.form.get('method','swanepoel')
        if method=='swanepoel':
            if 'T' not in material: raise ValueError('Swanepoel requires transmission')
            return jsonify(swanepoel(material['wavelength_nm'],material['T'],request.form.get('substrate_n','1.5')))
        if method=='tauc':
            if 'absorbance' not in material and 'T' not in material: raise ValueError('Tauc requires absorbance or transmission')
            return jsonify(tauc(material['wavelength_nm'],material.get('absorbance'),material.get('T'),request.form.get('thickness_nm','1000'),request.form.get('exponent','2'),(request.form.get('fit_lo_ev','2'),request.form.get('fit_hi_ev','3'))))
        raise ValueError('Unknown classical analysis method')
    except (ValueError,TypeError,KeyError) as exc: return jsonify({'error':str(exc)}),400

@app.post('/api/stack-simulate')
def stack_simulate():
    try:
        payload=request.get_json(force=True)
        w=np.asarray(payload['wavelength_nm'],dtype=float)
        if w.size>3000 or w.size<2: raise ValueError('Supply 2–3000 wavelengths')
        layers=payload.get('layers',[])
        if len(layers)>20: raise ValueError('Maximum 20 coherent layers')
        args=(w,layers,payload.get('ambient_n',1),payload.get('substrate_n',1.5),payload.get('angle_deg',0))
        result=stack_spectrum(*args,payload.get('polarization','s'))
        if payload.get('ellipsometry',False):
            result['ellipsometry']=ellipsometry(*args)
        def convert(v):
            if isinstance(v,np.ndarray):return v.tolist()
            if isinstance(v,dict):return {k:convert(x) for k,x in v.items()}
            return v
        return jsonify(convert(result))
    except (ValueError,TypeError,KeyError,OverflowError) as exc:return jsonify({'error':str(exc)}),400

@app.get('/api/reference')
def get_ref():
    try:
        s,(w,n,k)=samples()
    except (OSError, ValueError, KeyError) as exc:
        return jsonify({'error': 'Bundled R6G example unavailable: '+str(exc)}),503
    return jsonify({'wavelength_nm':w.tolist(),'n':n.tolist(),'k':k.tolist(),'samples':[{'name':v.name,'thickness_nm':v.thickness_nm,'substrate_n':v.substrate_n,'R':v.R.tolist(),'T':v.T.tolist()} for v in s]})


@app.post('/api/custom-fit')
def custom_fit():
    try:
        uploaded=request.files.get('spectrum')
        if not uploaded:raise ValueError('Choose a CSV, TSV or JSON spectrum file')
        structure=request.form.get('structure','supported')
        if structure not in ('supported','freestanding'): raise ValueError('This structure is not yet supported by the R/T inverse model')
        ns='1' if structure=='freestanding' else request.form.get('substrate_n','1')
        s=parse_spectrum(uploaded.filename,uploaded.read(),request.form.get('thickness_nm','1000'),ns,request.form.get('wavelength_unit','auto'),request.form.get('signal_unit','auto'))
        method=request.form.get('optimizer','least_squares')
        if method not in ('least_squares','multistart','differential_evolution','hho','hho_hybrid'):raise ValueError('Unknown optimizer')
        basis=int(request.form.get('basis','8'));budget=int(request.form.get('max_nfev','100'))
        if not 4<=basis<=25 or not 5<=budget<=500:raise ValueError('Invalid basis or budget')
        settings=FitSettings(n_basis=basis,optimizer=method,max_nfev=budget,fit_thickness=request.form.get('fit_thickness')=='true')
        initialization={'method':'default','warning':None}
        if request.form.get('initialize_from_swanepoel')=='true':
            if s.T is None:
                initialization={'method':'default','warning':'Swanepoel initialization requires transmission.'}
            elif float(ns)<=1:
                initialization={'method':'default','warning':'Swanepoel initialization requires a transparent substrate with n>1.'}
            else:
                try:
                    seed=swanepoel(s.wavelength_nm,s.T,float(ns))
                    settings.initial_n=(seed['wavelength_nm'],seed['n'])
                    if seed['thickness_nm'] is not None:
                        settings.initial_thickness_nm=seed['thickness_nm']
                        settings.fit_thickness=True
                    initialization={'method':'swanepoel','thickness_nm':seed['thickness_nm'],'warnings':seed['warnings']}
                except (ValueError,TypeError,IndexError) as exc:
                    initialization={'method':'default','warning':'Swanepoel unavailable: '+str(exc)}

        result=fit([s],settings)
        nf,kf=result['optical_constants'](s.wavelength_nm)
        pred=compare_prediction(s,nf,kf,result['thicknesses'][0])
        return jsonify({'wavelength_nm':s.wavelength_nm.tolist(),'n':nf.tolist(),'k':kf.tolist(),'R_measured':None if s.R is None else s.R.tolist(),'T_measured':None if s.T is None else s.T.tolist(),'R_fit':pred['R'].tolist(),'T_fit':pred['T'].tolist(),'thickness_nm':float(result['thicknesses'][0]),'rmse':float(pred['rmse']),'success':bool(result['result'].success),'message':str(result['result'].message),'initialization':initialization,'note':'User-supplied spectrum: no ground-truth n/k available. Single-spectrum inversion may be non-unique.'})
    except Exception as exc:return jsonify({'error':str(exc)}),400

@app.post('/api/fit')
def fit_api():
    try:
        req=request.get_json(silent=True)
        if not isinstance(req,dict):raise ValueError('JSON object required')
        s,(w,n_true,k_true)=samples()
        chosen=set(req.get('samples',[0,1,2,3]));s=[v for i,v in enumerate(s) if i in chosen]
        if not s:raise ValueError('Select at least one sample')
        method=req.get('optimizer','least_squares')
        if method not in ['least_squares','multistart','differential_evolution','hho','hho_hybrid']:raise ValueError('Unsupported optimizer')
        settings=FitSettings(n_basis=int(req.get('basis',10)),smoothness=float(req.get('smoothness',.1)),optimizer=method,fit_thickness=bool(req.get('fit_thickness',False)),max_nfev=int(req.get('max_nfev',100)),iterations=int(req.get('iterations',12)),population=int(req.get('population',12)))
        if not 4<=settings.n_basis<=25 or not 0<=settings.smoothness<=100 or not 5<=settings.max_nfev<=500 or not 1<=settings.iterations<=100 or not 5<=settings.population<=100:raise ValueError('Parameter outside allowed range')
        result=fit(s,settings);nf,kf=result['optical_constants'](w)
        predictions=[]
        for i,v in enumerate(s):
            p=compare_prediction(v,nf,kf,result['thicknesses'][i]);predictions.append({'name':v.name,'rmse':p['rmse'],'R':p['R'].tolist(),'T':p['T'].tolist(),'R_measured':v.R.tolist(),'T_measured':v.T.tolist(),'thickness_nm':float(result['thicknesses'][i])})
        return jsonify({'wavelength_nm':w.tolist(),'n':nf.tolist(),'k':kf.tolist(),'n_true':n_true.tolist(),'k_true':k_true.tolist(),'rmse_n':float(np.sqrt(np.mean((nf-n_true)**2))),'rmse_k':float(np.sqrt(np.mean((kf-k_true)**2))),'weighted_rmse':result['weighted_rmse'],'predictions':predictions,'message':str(result['result'].message),'success':bool(result['result'].success),'evaluations':int(result['result'].nfev)})
    except Exception as e:return jsonify({'error':str(e)}),400

if __name__=='__main__':
    app.run(host=os.environ.get('HOST','127.0.0.1'),port=int(os.environ.get('PORT','8765')),debug=False)
