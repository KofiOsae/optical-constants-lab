"""Multilayer inverse API. Kept separate from established single-film routes."""
import json
import numpy as np
from flask import Blueprint, request, jsonify
from custom_spectra import parse_spectrum
from multilayer_inverse import StackSpectrum, predict_stack
from multilayer_dispersion import fit_dispersive_stacks

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

@multilayer_api.post('/api/multilayer/joint-fit')
def multilayer_joint_fit():
    try:
        uploads=request.files.getlist('spectra')
        configs=json.loads(request.form.get('configurations','[]'))
        if not 1<=len(uploads)<=12 or not isinstance(configs,list) or len(configs)!=len(uploads):
            raise ValueError('Supply 1–12 spectra and matching stack configurations')
        samples=[]
        for upload,c in zip(uploads,configs):
            if not isinstance(c,dict):raise ValueError('Invalid configuration')
            d=float(c['thickness_nm'])
            spec=parse_spectrum(upload.filename,upload.read(),d,c.get('substrate_n',1.5),
                                c.get('wavelength_unit','auto'),c.get('signal_unit','auto'))
            angle=float(c.get('angle_deg',0.))
            pol=str(c.get('polarization','s')).lower()
            if not 0<=angle<85 or pol not in ('s','p'):
                raise ValueError('Angle or polarization invalid')
            w=spec.wavelength_nm
            sample=StackSpectrum(w,spec.R,spec.T,
                before=parse_layers(c.get('before',[]),w),
                after=parse_layers(c.get('after',[]),w),
                substrate_n=float(c.get('substrate_n',1.5)),
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
        out=fit_dispersive_stacks(samples,model=model,
            n_basis=int(request.form.get('n_basis',8)),
            smoothness=float(request.form.get('smoothness',.01)),
            oscillators=int(request.form.get('oscillators',1)),
            fit_thickness=fit_thickness,thickness_bounds=bounds,
            max_nfev=min(1000,max(10,int(request.form.get('max_nfev',250)))))
        grid=np.linspace(min(s.wavelength_nm.min() for s,d in samples),
                         max(s.wavelength_nm.max() for s,d in samples),201)
        n,k=out['optical_constants'](grid)
        predictions=[]
        for i,(s,d) in enumerate(samples):
            ni,ki=out['optical_constants'](s.wavelength_nm)
            r,t=predict_stack(s,ni,ki,out['thicknesses'][i])
            errors=[]
            if s.R is not None:errors.extend((r-s.R).tolist())
            if s.T is not None:errors.extend((t-s.T).tolist())
            predictions.append({'sample':i,'thickness_nm':float(out['thicknesses'][i]),
                'rmse':float(np.sqrt(np.mean(np.square(errors)))),
                'wavelength_nm':s.wavelength_nm.tolist(),
                'R_fit':r.tolist(),'T_fit':t.tolist()})
        return jsonify({'success':bool(out['result'].success),
            'message':str(out['result'].message),'model':model,
            'wavelength_nm':grid.tolist(),'n':n.tolist(),'k':k.tolist(),
            'weighted_rmse':out['weighted_rmse'],'predictions':predictions,
            'sensitivity':out['sensitivity'],
            'dispersion_parameters':out['dispersion_parameters'],
            'limitations':'Coherent isotropic layers; one unknown material; semi-infinite lossless substrate. No incoherent backside reflections.'})
    except (ValueError,TypeError,KeyError,IndexError,OverflowError) as exc:
        return jsonify({'error':str(exc)}),400
