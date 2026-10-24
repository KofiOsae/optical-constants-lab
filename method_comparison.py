from dataclasses import replace
from engine import fit,compare_prediction
from analysis_methods import swanepoel

def compare_methods(s, settings, substrate_n):
    outputs=[]
    seed=None
    try:
        if s.T is None: raise ValueError('Transmission required')
        seed=swanepoel(s.wavelength_nm,s.T,float(substrate_n))
        outputs.append(dict(method='Swanepoel',status='analytical',**seed))
    except (ValueError,TypeError,KeyError) as exc:
        outputs.append(dict(method='Swanepoel',status='unavailable',reason=str(exc)))
    for name,hybrid in [('B-spline + TMM',False),('Swanepoel + B-spline + TMM',True)]:
        if hybrid and seed is None:
            outputs.append(dict(method=name,status='unavailable',reason='No valid Swanepoel initialization'))
            continue
        cfg=replace(settings)
        cfg.initial_n=(seed['wavelength_nm'],seed['n']) if hybrid else None
        cfg.initial_thickness_nm=seed.get('thickness_nm') if hybrid and cfg.fit_thickness else None
        try:
            f=fit([s],cfg)
            n,k=f['optical_constants'](s.wavelength_nm)
            d=float(f['thicknesses'][0])
            p=compare_prediction(s,n,k,d)
            outputs.append(dict(method=name,status='converged' if f['result'].success else 'stopped',wavelength_nm=s.wavelength_nm.tolist(),n=n.tolist(),k=k.tolist(),thickness_nm=d,rmse=float(p['rmse']),R_fit=p['R'].tolist(),T_fit=p['T'].tolist(),evaluations=int(f['result'].nfev)))
        except (ValueError,TypeError) as exc:
            outputs.append(dict(method=name,status='failed',reason=str(exc)))
    return dict(methods=outputs,wavelength_nm=s.wavelength_nm.tolist(),R_measured=None if s.R is None else s.R.tolist(),T_measured=None if s.T is None else s.T.tolist())
