from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from flask import Flask, jsonify, request, send_from_directory, send_file
import os
import numpy as np
import pandas as pd
from engine import FitSettings,fit,fit_joint_two_stage,simulate_four,compare_prediction,fit_historical_logk,estimate_broadband_seed,phase2_hybrid_seed,fit_absorption_aware_blind
from phase1_fresnel import phase1_fit
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
from method_comparison import compare_methods

@app.get('/api/benchmark/scientific-acceptance')
def scientific_acceptance():
    """Independent reference check; archived reconstruction is not fitted here."""
    try:
        from historical_reference import archived_phase2
        from optical_inversion_package.inversion import coherent_film
        reference_data=archived_phase2()
        w=np.asarray(reference_data['wavelength_nm'],float)
        n=np.asarray(reference_data['n'],float)
        k=np.asarray(reference_data['k'],float)
        d=float(reference_data['thickness_nm'])
        raw=__import__('json').loads((ROOT/'data/original_plot_raw_7.json').read_text())
        series={v['label']:v for v in raw['data']['series']}
        R=np.asarray(series['R']['y'],float)
        T=np.asarray(series['T']['y'],float)
        R_fit,T_fit=coherent_film(w,n+1j*k,d)
        train=np.arange(len(w))%5!=0
        def rmse(a,b,mask):
            return float(np.sqrt(np.mean((a[mask]-b[mask])**2)))
        metrics={'R_train_rmse':rmse(R_fit,R,train),'T_train_rmse':rmse(T_fit,T,train),
                 'R_holdout_rmse':rmse(R_fit,R,~train),'T_holdout_rmse':rmse(T_fit,T,~train),
                 'R_all_rmse':rmse(R_fit,R,np.ones(len(w),bool)),
                 'T_all_rmse':rmse(T_fit,T,np.ones(len(w),bool))}
        # Compare the documented archival forward model with the production coherent-film model.
        model_R=np.asarray(reference_data['R_fit'],float)
        model_T=np.asarray(reference_data['T_fit'],float)
        metrics['forward_R_max_abs']=float(np.max(np.abs(R_fit-model_R)))
        metrics['forward_T_max_abs']=float(np.max(np.abs(T_fit-model_T)))
        # Reference n/k accuracy is a separate, explicitly identified comparison.
        ref=pd.read_csv(ROOT/'data/R6G_Singh.csv')
        nw=np.asarray(ref['wavelength_um'],float)*1000
        nref=np.interp(w,nw,np.asarray(ref['n'],float))
        kref=np.interp(w,nw,np.asarray(ref['k'],float))
        metrics['n_reference_mae']=float(np.mean(np.abs(n-nref)))
        metrics['k_reference_mae']=float(np.mean(np.abs(k-kref)))
        metrics['n_reference_rmse']=float(np.sqrt(np.mean((n-nref)**2)))
        metrics['k_reference_rmse']=float(np.sqrt(np.mean((k-kref)**2)))
        limits={'forward_R_max_abs':1e-8,'forward_T_max_abs':1e-8}
        passed=all(metrics[key]<=value for key,value in limits.items())
        return jsonify({'passed':passed,'metrics':metrics,'limits':limits,
                        'source':'Archived R6G free-standing-film Phase II reconstruction',
                        'reference_only':True,
                        'warning':'Forward-model parity is tested. This is not a blind n/k recovery accuracy test.'})
    except (ValueError,TypeError,KeyError,OSError) as exc:
        return jsonify({'error':str(exc)}),500

@app.get('/api/benchmark/phase2-archive')
def phase2_archive():
    from historical_reference import archived_phase2
    return jsonify(archived_phase2())

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
def home():return send_from_directory(ROOT/'web','workspace.html')

@app.get('/legacy')
def legacy():return send_from_directory(ROOT/'web','index.html')

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


@app.post('/api/compare-methods')
def compare_methods_api():
    try:
        uploaded=request.files.get('spectrum')
        if not uploaded: raise ValueError('Choose a measured R/T file')
        structure=request.form.get('structure','supported')
        if structure!='supported': raise ValueError('Swanepoel comparison requires film on a transparent substrate')
        ns=request.form.get('substrate_n','1.5')
        spec=parse_spectrum(uploaded.filename,uploaded.read(),request.form.get('thickness_nm','1000'),ns,request.form.get('wavelength_unit','auto'),request.form.get('signal_unit','auto'))
        lo=float(request.form.get('wavelength_min_nm',spec.wavelength_nm[0]))
        hi=float(request.form.get('wavelength_max_nm',spec.wavelength_nm[-1]))
        mask=(spec.wavelength_nm>=lo)&(spec.wavelength_nm<=hi)
        if mask.sum()<12: raise ValueError('At least 12 points required in selected range')
        spec.wavelength_nm=spec.wavelength_nm[mask]
        if spec.R is not None:spec.R=spec.R[mask]
        if spec.T is not None:spec.T=spec.T[mask]
        basis=int(request.form.get('basis','8'))
        if not 4<=basis<=25: raise ValueError('Basis count must be 4–25')
        cfg=FitSettings(n_basis=basis,optimizer=request.form.get('optimizer','least_squares'),max_nfev=min(500,max(5,int(request.form.get('max_nfev','100')))),fit_thickness=request.form.get('fit_thickness')=='true')
        return jsonify(compare_methods(spec,cfg,float(ns)))
    except (ValueError,TypeError,KeyError) as exc:return jsonify({'error':str(exc)}),400

@app.post('/api/joint-fit')
def joint_fit_api():
    """Fit 1-12 thin-film spectra; one-film fits require a transparent dispersion model."""
    try:
        import json
        uploads=request.files.getlist('spectra')
        configs=json.loads(request.form.get('configurations','[]'))
        if not 1<=len(uploads)<=12 or not isinstance(configs,list) or len(configs)!=len(uploads):
            raise ValueError('Supply 1-12 spectra and matching configurations')
        spectra=[]
        for upload,c in zip(uploads,configs):
            if not isinstance(c,dict):
                raise ValueError('Each configuration must be an object')
            spec=parse_spectrum(upload.filename,upload.read(),c.get('thickness_nm',1000),
                c.get('substrate_n',1.5),c.get('wavelength_unit','auto'),c.get('signal_unit','auto'))
            spec.angle_deg=float(c.get('angle_deg',0))
            spec.polarization=str(c.get('polarization','s')).lower()
            spec.sigma_R=float(c.get('sigma_R',.005))
            spec.sigma_T=float(c.get('sigma_T',.005))
            if not np.isfinite(spec.angle_deg) or not 0<=spec.angle_deg<85:
                raise ValueError('Angle must be in [0,85) degrees')
            if spec.polarization not in ('s','p'):
                raise ValueError('Polarization must be s or p')
            spec.name=str(c.get('name',upload.filename or 'spectrum'))[:100]
            spectra.append(spec)
        basis=int(request.form.get('basis',10))
        budget=int(request.form.get('max_nfev',180))
        smooth=float(request.form.get('smoothness',.1))
        tolerance=float(request.form.get('thickness_tolerance',.15))
        optimizer=request.form.get('optimizer','multistart')
        if not 4<=basis<=40:
            raise ValueError('Spline basis must be between 4 and 40')
        if not 5<=budget<=500:
            raise ValueError('Evaluation budget must be between 5 and 500')
        if not 0<=smooth<=100:
            raise ValueError('Regularization must be between 0 and 100')
        if not 0<tolerance<=.5:
            raise ValueError('Thickness tolerance must be between 0 and 0.5')
        if optimizer not in ('least_squares','multistart','differential_evolution','hho','hho_hybrid'):
            raise ValueError('Unknown optimizer')
        n_min=float(request.form.get('n_min',0.05)); n_max=float(request.form.get('n_max',8.0)); k_max=float(request.form.get('k_max',6.0))
        if not (np.isfinite(n_min) and np.isfinite(n_max) and np.isfinite(k_max) and 0<n_min<n_max<=30 and 0<k_max<=30):
            raise ValueError('Require finite bounds 0 < n_min < n_max <= 30 and 0 < k_max <= 30')
        settings=FitSettings(n_basis=basis,max_nfev=budget,smoothness=smooth,
            optimizer=optimizer,fit_thickness=request.form.get('fit_thickness')=='true',
            thickness_tolerance=tolerance,n_min=n_min,n_max=n_max,k_max=k_max)
        method=request.form.get('joint_method','two_stage')
        if method not in ('two_stage','direct','lorentz','drude_lorentz','cauchy','sellmeier'): raise ValueError('Unknown joint method')
        if len(spectra)==1 and method not in ('cauchy','sellmeier'):
            raise ValueError('Single-film fitting currently supports Cauchy or Sellmeier only')
        holdout=int(request.form.get('holdout_index','-1'))
        if holdout < -1 or holdout >= len(spectra): raise ValueError('Invalid holdout index')
        if holdout >= 0 and len(spectra)<3: raise ValueError('Holdout validation needs at least 3 films')
        training=[sp for i,sp in enumerate(spectra) if i!=holdout]
        # Stable fingerprint of the numerical training problem, independent of filenames.
        import hashlib
        fingerprint=hashlib.sha256()
        for spec in training:
            for values in (spec.wavelength_nm,spec.R,spec.T):
                if values is None:
                    fingerprint.update(b'NONE')
                else:
                    arr=np.asarray(values,dtype='<f8')
                    fingerprint.update(str(arr.shape).encode('ascii'))
                    fingerprint.update(arr.tobytes())
            fingerprint.update(np.asarray([spec.thickness_nm,spec.substrate_n,spec.ambient_n,spec.angle_deg,spec.sigma_R,spec.sigma_T],dtype='<f8').tobytes())
            fingerprint.update(str(spec.polarization).encode('utf-8'))
        training_fingerprint=fingerprint.hexdigest()

        if method in ('lorentz','drude_lorentz'):
            from scientific_diagnostics import fit_causal
            outcome=fit_causal(training,method,int(request.form.get('oscillators',1)),budget,settings.fit_thickness,multistart=int(request.form.get('causal_starts','1')),warm_start=(json.loads(request.form['causal_warm_start']) if request.form.get('causal_warm_start') else None))
        elif method in ('cauchy','sellmeier'):
            from dispersion_models import fit_transparent_spectra
            outcome=fit_transparent_spectra(training,method,budget,settings.fit_thickness)
        else:
            outcome=fit_joint_two_stage(training,settings) if method=='two_stage' else fit(training,settings)
        grid=np.unique(np.concatenate([s.wavelength_nm for s in training]))
        n,k=outcome['optical_constants'](grid)
        predictions=[]
        for i,spec in enumerate(spectra):
            ni,ki=outcome['optical_constants'](spec.wavelength_nm)
            fitted_index=i if holdout<0 or i<holdout else i-1
            thickness=spec.thickness_nm if i==holdout else outcome['thicknesses'][fitted_index]
            p=compare_prediction(spec,ni,ki,thickness)
            predictions.append(dict(name=spec.name,wavelength_nm=spec.wavelength_nm.tolist(),
                R_measured=None if spec.R is None else spec.R.tolist(),
                T_measured=None if spec.T is None else spec.T.tolist(),
                R_fit=p['R'].tolist(),T_fit=p['T'].tolist(),rmse=p['rmse'],
                thickness_nm=float(thickness),held_out=bool(i==holdout),
                angle_deg=spec.angle_deg,polarization=spec.polarization))
        from scientific_diagnostics import jacobian_report
        diagnostics=jacobian_report(outcome['result'])
        return jsonify(dict(wavelength_nm=grid.tolist(),n=n.tolist(),k=k.tolist(),
            predictions=predictions,weighted_rmse=outcome['weighted_rmse'],training_fingerprint=training_fingerprint,stages=outcome.get('stages'),
            sensitivity=diagnostics,holdout_index=holdout,
            success=bool(outcome['result'].success),message=str(outcome['result'].message),
            note='Shared n,k; coherent uniform films on semi-infinite lossless substrates.'))
    except (ValueError,TypeError,KeyError,OverflowError) as exc:
        return jsonify({'error':str(exc)}),400

@app.post('/api/custom-fit')
def custom_fit():
    try:
        uploaded=request.files.get('spectrum')
        if not uploaded:raise ValueError('Choose a CSV, TSV or JSON spectrum file')
        structure=request.form.get('structure','supported')
        if structure not in ('supported','freestanding'): raise ValueError('This structure is not yet supported by the R/T inverse model')
        ns='1' if structure=='freestanding' else request.form.get('substrate_n','1')
        s=parse_spectrum(uploaded.filename,uploaded.read(),request.form.get('thickness_nm','1000'),ns,request.form.get('wavelength_unit','auto'),request.form.get('signal_unit','auto'))
        # Restrict the numerical inversion itself to the user-selected wavelengths.
        range_lo=request.form.get('wavelength_min_nm')
        range_hi=request.form.get('wavelength_max_nm')
        if (range_lo is not None or range_hi is not None) and request.form.get('reconstruction_model')!='historical_logk':
            lo=float(range_lo) if range_lo not in (None,'') else float(s.wavelength_nm[0])
            hi=float(range_hi) if range_hi not in (None,'') else float(s.wavelength_nm[-1])
            if not np.isfinite(lo) or not np.isfinite(hi) or lo>=hi:
                raise ValueError('Invalid wavelength range: minimum must be below maximum')
            mask=(s.wavelength_nm>=lo)&(s.wavelength_nm<=hi)
            if np.count_nonzero(mask)<max(8,int(request.form.get('basis','8'))+2):
                raise ValueError('Selected wavelength range contains too few points for the chosen spline basis')
            s.wavelength_nm=s.wavelength_nm[mask]
            if s.R is not None:s.R=s.R[mask]
            if s.T is not None:s.T=s.T[mask]
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

        if request.form.get('reconstruction_model')=='phase1':
            if s.R is None or s.T is None:
                raise ValueError('Phase I requires both R and T')
            result=phase1_fit(s.wavelength_nm,s.R,s.T,max_nfev=budget,substrate_n=float(ns),thickness_estimate_nm=s.thickness_nm)
            result['initialization']={'method':'documented_phase1_multistart'}
            result['note']='Phase I weak-absorption Fresnel model; supports air/film/air and film on semi-infinite transparent substrate. Original R6G benchmark is freestanding.'
            return jsonify(result)
        if request.form.get('reconstruction_model')=='historical_logk' and (s.R is None or s.T is None):
            raise ValueError('Phase II requires both reflectance and transmittance')
        if request.form.get('reconstruction_model')=='historical_logk':
            phase1_lo=request.form.get('wavelength_min_nm')
            phase1_hi=request.form.get('wavelength_max_nm')
            if phase1_lo in (None,'') or phase1_hi in (None,''):
                raise ValueError('Select Phase I initialization window using the wavelength sliders')
            phase1_lo=float(phase1_lo);phase1_hi=float(phase1_hi)
            if not np.isfinite(phase1_lo) or not np.isfinite(phase1_hi) or phase1_lo>=phase1_hi:
                raise ValueError('Invalid Phase I initialization range')
            mask=(s.wavelength_nm>=phase1_lo)&(s.wavelength_nm<=phase1_hi)
            phase1=None
            approved_text=request.form.get('approved_phase1','').strip()
            if approved_text:
                import json
                approved=json.loads(approved_text)
                p=approved.get('parameters')
                if not isinstance(p,dict) or not all(key in p for key in ('a','b','c','q0','q1','d_nm')):
                    raise ValueError('Approved Phase I must contain all six fitted parameters')
                if approved.get('R_measured') != s.R[mask].tolist() or approved.get('T_measured') != s.T[mask].tolist():
                    raise ValueError('Phase I input data changed; rerun Phase I')
                if approved.get('wavelength_nm') != s.wavelength_nm[mask].tolist():
                    raise ValueError('Phase I wavelength selection changed; rerun Phase I')
                if approved.get('substrate_n') != float(ns) or approved.get('structure') != structure:
                    raise ValueError('Sample geometry changed; rerun Phase I')
                p={key:float(p[key]) for key in ('a','b','c','q0','q1','d_nm')}
                if not all(np.isfinite(v) for v in p.values()) or p['d_nm']<=0:
                    raise ValueError('Invalid approved Phase I parameters')
                phase1={'parameters':p,'thickness_nm':p['d_nm'],'wavelength_nm':s.wavelength_nm[mask].tolist(),'approved':True}
            if np.count_nonzero(mask)<8:
                raise ValueError('Phase I initialization range requires at least 8 measured R/T points')
            if np.count_nonzero(mask)>=8:
                if phase1 is None:
                    phase1=phase1_fit(s.wavelength_nm[mask],s.R[mask],s.T[mask],max_nfev=budget,substrate_n=float(ns),thickness_estimate_nm=s.thickness_nm)
                s.thickness_nm=phase1['thickness_nm']
                p=phase1['parameters']
                x=800.0/s.wavelength_nm
                n_seed=p['a']+p['b']*x*x+p['c']*x**4
                k_seed=np.exp(np.clip(p['q0']+p['q1']*(x-1),-24,2))
                broadband=estimate_broadband_seed(s,n_seed,k_seed,max_nfev=max(20,min(70,budget)))
                hybrid=phase2_hybrid_seed(s.wavelength_nm,s.wavelength_nm,broadband['n'],broadband['k'],p)
                historical=fit_historical_logk(s,max_nfev=budget,warm_start=hybrid)
            return jsonify({'wavelength_nm':s.wavelength_nm.tolist(),
                'n':historical['n'].tolist(),'k':historical['k'].tolist(),
                'R_measured':s.R.tolist(),'T_measured':s.T.tolist(),
                'R_fit':historical['R'].tolist(),'T_fit':historical['T'].tolist(),
                'thickness_nm':float(s.thickness_nm),
                'rmse':historical['train_rmse'],
                'train_rmse':historical['train_rmse'],'holdout_rmse':historical['holdout_rmse'],
                'split':['train' if v else 'holdout' for v in historical['train_mask']],
                'success':bool(historical['result'].success),
                'message':str(historical['result'].message),'phase2_nfev':historical['nfev'],'phase2_candidate_count':historical['candidate_count'],'phase2_regularized_objective':historical['regularized_objective'],
                'initialization':{'method':'sample_specific_broadband_plus_phase1_hybrid','broadband_nfev':broadband['nfev'],'broadband_success':broadband['success'],'transition_nm':690.0},'phase1':phase1,
                'note':'Sample-specific broadband R/T initializer plus Phase I above 690 nm; Phase II holds out every fifth wavelength. Historical R6G archive is not used for user samples.'})
        if request.form.get('reconstruction_model')=='blind_absorption_aware':
            if s.R is None or s.T is None:
                raise ValueError('Blind absorption-aware inversion requires both R and T')
            if settings.fit_thickness:
                raise ValueError('Blind initializer requires fixed thickness; use a thickness sweep for sensitivity')
            if request.form.get('initialize_from_swanepoel')=='true':
                raise ValueError('Select either Swanepoel or blind absorption-aware initialization')
            threshold_text=request.form.get('transparent_min_nm','')
            threshold=float(threshold_text) if threshold_text.strip() else None
            if threshold is not None and not np.isfinite(threshold):
                raise ValueError('Invalid transparent-window threshold')
            settings.smoothness=float(request.form.get('smoothness','0.01'))
            if not 0<=settings.smoothness<=100:
                raise ValueError('Smoothness must be between 0 and 100')
            result=fit_absorption_aware_blind(s,settings,transparent_min_nm=threshold)
            return jsonify({'wavelength_nm':s.wavelength_nm.tolist(),
                'n':result['n'].tolist(),'k':result['k'].tolist(),
                'R_measured':s.R.tolist(),'T_measured':s.T.tolist(),
                'R_fit':result['R_fit'].tolist(),'T_fit':result['T_fit'].tolist(),
                'thickness_nm':float(result['thickness_nm']),
                'rmse':float(np.sqrt((result['R_rmse']**2+result['T_rmse']**2)/2)),
                'R_rmse':result['R_rmse'],'T_rmse':result['T_rmse'],
                'success':result['success'],'message':'Blind staged inversion completed',
                'initialization':{'method':'absorption_aware_blind','starts':result['candidates'],
                    'transparent_threshold_nm':result['transparent_threshold_nm'],'nfev':result['nfev']},
                'note':'Data-only initializer; thickness fixed. Optical constants may be nonunique.'})
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

@app.post('/api/export-causal-workbook')
def export_causal_workbook():
    try:
        from model_workbook import build_model_workbook
        payload=request.get_json(force=True)
        model=payload['model']
        coefficients=payload['coefficients']
        wavelengths=payload['wavelength_nm']
        if len(wavelengths)>10000 or len(coefficients)>12:
            raise ValueError('Workbook size exceeds allowed limits')
        import numpy as np
        if not np.all(np.isfinite(np.asarray(coefficients,dtype=float))) or not np.all(np.isfinite(np.asarray(wavelengths,dtype=float))) or np.any(np.asarray(wavelengths,dtype=float)<=0):
            raise ValueError('Wavelengths and coefficients must be finite and wavelengths positive')
        book=build_model_workbook(model,coefficients,wavelengths,payload.get('metadata',{}))
        return send_file(book,as_attachment=True,download_name='fitted_optical_model.xlsx',mimetype='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')
    except (ValueError,TypeError,KeyError) as exc:
        return jsonify({'error':str(exc)}),400
