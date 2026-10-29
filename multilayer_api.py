"""Multilayer inverse API. Kept separate from established single-film routes."""
import json
import numpy as np
from flask import Blueprint, request, jsonify
from custom_spectra import parse_spectrum, apply_rt_noise
from material_import import inspect_material
from multilayer_inverse import StackSpectrum, predict_stack
from multilayer_dispersion import fit_dispersive_stacks
from kk_sskk import kk_finite_band, sskk_finite_band

multilayer_api=Blueprint('multilayer_api',__name__)

def parse_layers(layers,w):
    if not isinstance(layers,list) or len(layers)>12:
        raise ValueError('Provide up to 12 known layers')
    result=[]
    for item in layers:
        if not isinstance(item,dict):raise ValueError('Invalid layer')
        d=float(item['thickness_nm'])
        if not np.isfinite(d) or d<0:raise ValueError('Invalid layer thickness')
        if 'wavelength_nm' in item:
            grid=np.asarray(item['wavelength_nm'],float)
            n=np.asarray(item['n'],float)
            k=np.asarray(item.get('k',np.zeros_like(n)),float)
            if len(grid)<2 or np.any(np.diff(grid)<=0) or grid[0]>w.min() or grid[-1]<w.max():
                raise ValueError('Known-layer RI table must cover measurement spectrum')
            if n.shape!=grid.shape or k.shape!=grid.shape:
                raise ValueError('Known-layer RI array length mismatch')
            value=np.interp(w,grid,n)+1j*np.interp(w,grid,k)
        else:
            value=complex(float(item['n']),float(item.get('k',0.)))
        result.append((value,d))
    return tuple(result)

def substrate_index(configuration,w):
    table=configuration.get('substrate_ri')
    if table is None:
        return complex(float(configuration.get('substrate_n',1.5)),float(configuration.get('substrate_k',0)))
    if not isinstance(table,dict):
        raise ValueError('Substrate RI table must be an object')
    return parse_layers([dict(table,thickness_nm=0.)],w)[0][0]

@multilayer_api.post('/api/multilayer/preview')
def multilayer_preview():
    try:
        upload=request.files.get('spectrum')
        if upload is None:
            raise ValueError('Upload a measured R/T spectrum')
        data=inspect_material(upload.filename,upload.read(),
            request.form.get('wavelength_unit','auto'),
            request.form.get('signal_unit','auto'))
        columns=data['columns']
        if 'R' not in columns or 'T' not in columns:
            raise ValueError('Both R and T are required')
        if data['point_count']<8:
            raise ValueError('At least 8 wavelength points required')
        return jsonify({'wavelength_nm':data['wavelength_nm'],
            'R':columns['R'],'T':columns['T'],
            'point_count':data['point_count'],
            'wavelength_range_nm':data['wavelength_range_nm']})
    except (ValueError,TypeError,KeyError) as exc:
        return jsonify({'error':str(exc)}),400

@multilayer_api.post('/api/multilayer/joint-fit')
def multilayer_joint_fit():
    try:
        uploads=request.files.getlist('spectra')
        configs=json.loads(request.form.get('configurations','[]'))
        if not 1<=len(uploads)<=12 or not isinstance(configs,list) or len(configs)!=len(uploads):
            raise ValueError('Supply 1–12 spectra and matching stack configurations')
        samples=[]
        sigma=float(request.form.get('noise_sigma','0') or '0')
        seed=int(request.form.get('noise_seed','2026') or '2026')
        for i,(upload,c) in enumerate(zip(uploads,configs)):
            if not isinstance(c,dict):raise ValueError('Invalid configuration')
            d=float(c['thickness_nm'])
            spec=parse_spectrum(upload.filename,upload.read(),d,c.get('substrate_n',1.5),
                                c.get('wavelength_unit','auto'),c.get('signal_unit','auto'))
            limits=c.get('fit_range_nm')
            if limits is not None:
                lo,hi=map(float,limits)
                if not np.isfinite(lo) or not np.isfinite(hi) or lo>=hi:
                    raise ValueError('Invalid fitting wavelength range')
                mask=(spec.wavelength_nm>=lo)&(spec.wavelength_nm<=hi)
                if np.count_nonzero(mask)<8:
                    raise ValueError('Fitting wavelength range must contain at least 8 measured points')
                spec.wavelength_nm=spec.wavelength_nm[mask]
                if spec.R is not None: spec.R=spec.R[mask]
                if spec.T is not None: spec.T=spec.T[mask]
            spec=apply_rt_noise(spec,sigma,seed,i)
            angle=float(c.get('angle_deg',0.))
            pol=str(c.get('polarization','s')).lower()
            if not 0<=angle<85 or pol not in ('s','p'):
                raise ValueError('Angle or polarization invalid')
            w=spec.wavelength_nm
            sample=StackSpectrum(w,spec.R,spec.T,
                roughness_nm=tuple(float(v) for v in c['interface_mixing_nm']) if c.get('interface_mixing_nm') is not None else None,
                before=parse_layers(c.get('before',[]),w),
                after=parse_layers(c.get('after',[]),w),
                substrate_n=substrate_index(c,w),
                substrate_thickness_nm=float(c['substrate_thickness_nm']) if c.get('substrate_thickness_nm') not in (None,'') else None,
                backing_n=float(c.get('backing_n',1.)),
                ambient_n=float(c.get('ambient_n',1.)),
                angle_deg=angle,polarization=pol,
                sigma_R=float(c.get('sigma_R',.005)),
                sigma_T=float(c.get('sigma_T',.005)))
            samples.append((sample,d))
        model=str(request.form.get('model','spline'))
        fit_thickness=request.form.get('fit_thickness','false').lower()=='true'
        bounds=None
        if fit_thickness:
            bounds=[tuple(map(float,c['thickness_bounds_nm'])) for c in configs]
        initial=None
        warm_raw=request.form.get('warm_start','').strip()
        fit_low=min(float(s.wavelength_nm.min()) for s,d in samples)
        fit_high=max(float(s.wavelength_nm.max()) for s,d in samples)
        basis=int(request.form.get('n_basis',8))
        oscillator_count=int(request.form.get('oscillators',1))
        expanded_warm=None
        expanded_raw=request.form.get('warm_start_spectrum','').strip()
        if expanded_raw:
            expanded_warm=json.loads(expanded_raw)
            if not isinstance(expanded_warm,dict):
                raise ValueError('Stage-one warm transfer must contain a fitted optical spectrum')
            if expanded_warm.get('source')!='previous_multilayer_fit':
                raise ValueError('Expanded warm transfer must originate from a previous multilayer fit')
            if expanded_warm.get('model')!=model:
                raise ValueError('Stage-one and expanded fit must use the same model')
        if warm_raw and expanded_warm is not None:
            raise ValueError('Choose an exact warm start or an expanded spectrum transfer, not both')
        if warm_raw:
            warm=json.loads(warm_raw)
            if not isinstance(warm,dict) or warm.get('model')!=model:
                raise ValueError('Warm-start model must match selected dispersion model')
            if int(warm.get('n_basis',-1))!=basis or int(warm.get('oscillators',-1))!=oscillator_count:
                raise ValueError('Warm-start basis and oscillator settings must match')
            domain=warm.get('wavelength_range_nm',[])
            if len(domain)!=2 or not np.allclose(domain,[fit_low,fit_high],rtol=0,atol=1e-5):
                raise ValueError('Warm-start wavelength domain must match the current fitting domain')
            initial=warm.get('dispersion_parameters')
        out=fit_dispersive_stacks(samples,model=model,
            n_basis=basis,
            smoothness=float(request.form.get('smoothness',.01)),
            oscillators=oscillator_count,
            fit_thickness=fit_thickness,thickness_bounds=bounds,
            optical_bounds=tuple(float(request.form.get(key,default)) for key,default in [('n_min',.5),('n_max',6.),('k_min',0.),('k_max',5.)]),
            max_nfev=min(1000,max(10,int(request.form.get('max_nfev',250)))),
            initial_parameters=initial,initial_reference=expanded_warm,
            parameter_controls=json.loads(request.form['oscillator_parameters']) if request.form.get('oscillator_parameters') else None,
            fit_interface_mixing=request.form.get('fit_interface_mixing','false')=='true',
            interface_mixing_bounds=(float(request.form.get('mixing_min_nm',0)),float(request.form.get('mixing_max_nm',10))))
        grid=np.linspace(min(s.wavelength_nm.min() for s,d in samples),
                         max(s.wavelength_nm.max() for s,d in samples),201)
        n,k=out['optical_constants'](grid)
        diagnostic=None
        mode=request.form.get('dispersion_diagnostic','none')
        if mode not in ('none','kk','sskk'):raise ValueError('Unsupported dispersion diagnostic')
        if mode!='none':
            eps1=n*n-k*k
            eps2=2*n*k
            if mode=='kk':
                diagnostic=kk_finite_band(grid,eps2,epsilon_inf=float(request.form.get('epsilon_inf',1.)))
            else:
                aw=float(request.form['anchor_wavelength_nm'])
                ae=float(request.form['anchor_epsilon1'])
                diagnostic=sskk_finite_band(grid,eps2,aw,ae)
            diagnostic['mode']=mode
            diagnostic['epsilon1_fitted']=np.interp(diagnostic['wavelength_nm'],grid,eps1).tolist()
            diagnostic['status']='finite_band_diagnostic_only'
        predictions=[]
        for i,(s,d) in enumerate(samples):
            ni,ki=out['optical_constants'](s.wavelength_nm)
            if out.get('fitted_interface_mixing_nm') is not None:
                import dataclasses
                count=len(s.before)+len(s.after)+2+(1 if s.substrate_thickness_nm is not None else 0)
                s=dataclasses.replace(s,roughness_nm=(out['fitted_interface_mixing_nm'],)*count)
            r,t=predict_stack(s,ni,ki,out['thicknesses'][i])
            errors=[]
            if s.R is not None:errors.extend((r-s.R).tolist())
            if s.T is not None:errors.extend((t-s.T).tolist())
            predictions.append({'sample':i,'thickness_nm':float(out['thicknesses'][i]),
                'nominal_thickness_nm':float(d),
                'thickness_fitted':bool(fit_thickness),
                'thickness_bounds_nm':list(map(float,configs[i]['thickness_bounds_nm'])) if fit_thickness else None,
                'rmse':float(np.sqrt(np.mean(np.square(errors)))),
                'wavelength_nm':s.wavelength_nm.tolist(),
                'R_fit':r.tolist(),'T_fit':t.tolist(),
                'R_measured':s.R.tolist() if s.R is not None else None,
                'T_measured':s.T.tolist() if s.T is not None else None})
        return jsonify({'success':bool(out['result'].success),
            'message':str(out['result'].message),'model':model,
            'fit_thickness':bool(fit_thickness),
            'noise':{'sigma':sigma,'seed':seed,'applied':sigma>0,
                     'clipped_to_physical_range':sigma>0},
            'warm_start_used':bool(out['warm_start_used']),
            'warm_start_transferred':bool(out['warm_start_transferred']),
            'warm_start_package':{'model':model,'n_basis':basis,
                'oscillators':oscillator_count,'wavelength_range_nm':[fit_low,fit_high],
                'dispersion_parameters':out['dispersion_parameters']},
            'wavelength_nm':grid.tolist(),'n':n.tolist(),'k':k.tolist(),
            'weighted_rmse':out['weighted_rmse'],'fitted_interface_mixing_nm':out.get('fitted_interface_mixing_nm'),'predictions':predictions,
            'dispersion_diagnostic':diagnostic,
            'sensitivity':out['sensitivity'],
            'dispersion_parameters':out['dispersion_parameters'],
            'optical_bounds':out.get('optical_bounds'),
            'multistart':out.get('multistart'),
            'ambiguity':out.get('ambiguity'),
            'model_note':out.get('model_note',''),
            'limitations':'Coherent isotropic layers; one unknown material; semi-infinite lossless substrate. No incoherent backside reflections.'})
    except (ValueError,TypeError,KeyError,IndexError,OverflowError) as exc:
        return jsonify({'error':str(exc)}),400

@multilayer_api.post('/api/multilayer/forward-validation')
def multilayer_forward_validation():
    """Independent reference-RI forward prediction against uploaded R/T."""
    try:
        spectrum=request.files.get('spectrum')
        reference=request.files.get('reference')
        if spectrum is None or reference is None:
            raise ValueError('Upload both R/T spectrum and reference n,k CSV')
        cfg=json.loads(request.form.get('configuration','{}'))
        d=float(cfg['thickness_nm'])
        if not np.isfinite(d) or d<=0:raise ValueError('Thickness must be positive')
        spec=parse_spectrum(spectrum.filename,spectrum.read(),d,cfg.get('substrate_n',1.5),
                            cfg.get('wavelength_unit','auto'),cfg.get('signal_unit','auto'))
        import io,csv
        raw=reference.read().decode('utf-8-sig')
        reader=csv.DictReader(io.StringIO(raw),delimiter=';' if ';' in raw.splitlines()[0] else ',')
        if not reader.fieldnames:raise ValueError('Reference CSV has no header')
        names={h.strip().lower():h for h in reader.fieldnames}
        wkey=names.get('wavelength_nm') or names.get('wavelength_um')
        if not wkey or 'n' not in names or 'k' not in names:
            raise ValueError('Reference requires wavelength_nm or wavelength_um, n, k')
        rows=list(reader)
        if len(rows)<2:raise ValueError('Reference needs at least two rows')
        rw=np.array([float(row[wkey]) for row in rows])
        if wkey==names.get('wavelength_um'):rw*=1000
        rn=np.array([float(row[names['n']]) for row in rows])
        rk=np.array([float(row[names['k']]) for row in rows])
        if not np.all(np.isfinite(rw)) or not np.all(np.isfinite(rn)) or not np.all(np.isfinite(rk)) or np.any(np.diff(rw)<=0) or np.any(rn<=0) or np.any(rk<0):
            raise ValueError('Reference wavelength must increase strictly; n>0, k>=0 and all values finite')
        w=np.asarray(spec.wavelength_nm,float)
        if w.min()<rw.min()-1e-7 or w.max()>rw.max()+1e-7:
            raise ValueError('Reference RI must cover the full measured wavelength interval; extrapolation is disabled')
        ni=np.interp(w,rw,rn);ki=np.interp(w,rw,rk)
        sample=StackSpectrum(w,spec.R,spec.T,
            before=parse_layers(cfg.get('before',[]),w),
            after=parse_layers(cfg.get('after',[]),w),
            substrate_n=substrate_index(cfg,w),
            substrate_thickness_nm=float(cfg['substrate_thickness_nm']) if cfg.get('substrate_thickness_nm') not in (None,'') else None,
            backing_n=float(cfg.get('backing_n',1.)),
            ambient_n=float(cfg.get('ambient_n',1.)),
            angle_deg=float(cfg.get('angle_deg',0.)),
            polarization=str(cfg.get('polarization','s')).lower())
        r,t=predict_stack(sample,ni,ki,d)
        def metric(obs,pred):
            return None if obs is None else float(np.sqrt(np.mean((np.asarray(obs)-pred)**2)))
        errors=[(r-spec.R)] if spec.R is not None else []
        if spec.T is not None:errors.append(t-spec.T)
        return jsonify({'wavelength_nm':w.tolist(),'R_predicted':r.tolist(),'T_predicted':t.tolist(),
            'R_measured':None if spec.R is None else spec.R.tolist(),
            'T_measured':None if spec.T is None else spec.T.tolist(),
            'R_residual':None if spec.R is None else (r-spec.R).tolist(),
            'T_residual':None if spec.T is None else (t-spec.T).tolist(),
            'rmse_R':metric(spec.R,r),'rmse_T':metric(spec.T,t),
            'rmse_combined':float(np.sqrt(np.mean(np.concatenate(errors)**2))),
            'note':'No fitting. Reference interpolated only within measured range; coherent semi-infinite-substrate model.'})
    except (ValueError,TypeError,KeyError,IndexError,UnicodeError) as exc:
        return jsonify({'error':str(exc)}),400
