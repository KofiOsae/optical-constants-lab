"""Transparent-material dispersion models. Wavelength input in nanometres."""
import numpy as np

def cauchy(wavelength_nm, A, B, C=0.0):
    w=np.asarray(wavelength_nm,dtype=float)/1000
    if np.any(w<=0):raise ValueError('Wavelength must be positive')
    return A+B/w**2+C/w**4

def sellmeier(wavelength_nm, coefficients):
    w=np.asarray(wavelength_nm,dtype=float)/1000
    if np.any(w<=0):raise ValueError('Wavelength must be positive')
    if len(coefficients)!=6:raise ValueError('Sellmeier needs B1,C1,B2,C2,B3,C3')
    n2=np.ones_like(w)
    for B,C in zip(coefficients[::2],coefficients[1::2]):
        denominator=w*w-C
        if np.any(np.abs(denominator)<1e-12):
            raise ValueError('Sellmeier pole lies at a requested wavelength')
        n2+=B*w*w/denominator
    if np.any(n2<=0):raise ValueError('Sellmeier n squared is nonpositive in the requested range')
    return np.sqrt(n2)

def fit_cauchy(wavelength_nm, n, order=3):
    """Fit n=A+B/lambda_um^2+C/lambda_um^4 using linear least squares."""
    w=np.asarray(wavelength_nm,float)/1000
    n=np.asarray(n,float)
    if w.ndim!=1 or n.shape!=w.shape or len(w)<order or order not in (2,3):
        raise ValueError('Invalid input or Cauchy order')
    if not np.all(np.isfinite(w)) or not np.all(np.isfinite(n)) or np.any(w<=0):
        raise ValueError('Finite positive wavelengths and finite indices required')
    X=np.column_stack([w**(-2*j) for j in range(order)])
    coeff,_,rank,_=np.linalg.lstsq(X,n,rcond=None)
    if rank!=order:raise ValueError('Cauchy fit is rank deficient')
    return dict(coefficients=coeff.tolist(),rmse=float(np.sqrt(np.mean((X@coeff-n)**2))))

def fit_transparent_spectra(spectra, model='cauchy', max_nfev=180, fit_thickness=False):
    """Joint R/T inversion for transparent films; k=0. Sellmeier C values fixed to avoid unconstrained poles."""
    from scipy.optimize import least_squares
    from engine import coherent_film
    if model not in ('cauchy','sellmeier'):raise ValueError('Unknown transparent dispersion model')
    nom=np.array([s.thickness_nm for s in spectra],float)
    if model=='cauchy':
        initial=np.array([1.5,.005,.0001])
        lower=np.array([1.0,-.5,-.2])
        upper=np.array([4.0,.5,.2])
    else:
        # Fixed resonance wavelengths (UV/IR), fitting oscillator strengths only.
        initial=np.array([.7,.4,.9])
        lower=np.array([0.,0.,0.])
        upper=np.array([5.,5.,5.])
    poles=np.array([.0684043**2,.1162414**2,9.896161**2])
    def optical(w,p):
        if model=='cauchy':n=cauchy(w,*p[:3])
        else:
            params=[v for pair in zip(p[:3],poles) for v in pair]
            n=sellmeier(w,params)
        return n,np.zeros_like(n)
    def residual(p):
        output=[]
        for i,s in enumerate(spectra):
            n,k=optical(s.wavelength_nm,p)
            d=p[len(initial)+i] if fit_thickness else s.thickness_nm
            R,T=coherent_film(s.wavelength_nm,n+1j*k,d,s.substrate_n,s.ambient_n,s.angle_deg,s.polarization)
            if s.R is not None:output.extend((R-s.R)/s.sigma_R)
            if s.T is not None:output.extend((T-s.T)/s.sigma_T)
        return np.asarray(output)
    if fit_thickness:
        initial=np.r_[initial,nom]
        lower=np.r_[lower,nom*.95]
        upper=np.r_[upper,nom*1.05]
    fitted=least_squares(residual,initial,bounds=(lower,upper),max_nfev=max_nfev,x_scale='jac')
    p=fitted.x.copy()
    return dict(result=fitted,optical_constants=lambda w:optical(w,p),
        thicknesses=p[3:] if fit_thickness else nom,
        weighted_rmse=float(np.sqrt(np.mean(residual(p)**2))),
        stages={'selected':model,'coefficients':p[:3].tolist(),
                'sellmeier_fixed_C_um2':poles.tolist() if model=='sellmeier' else None,
                'assumption':'k=0 (transparent film); Sellmeier resonance positions fixed to fused-silica reference values' if model=='sellmeier' else 'k=0 (transparent film)'},
        model=model)
