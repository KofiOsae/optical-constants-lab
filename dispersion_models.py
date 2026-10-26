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
