"""Reference multi-spectrum optical-constant inversion engine (not integrated into Workbench).
Time dependence exp(-i omega t); n_complex=n+i*k; coherent single film,
semi-infinite ambient and substrate. R,T are power fractions, not percentages.
"""
from dataclasses import dataclass
import numpy as np
from scipy.interpolate import BSpline
from scipy.optimize import least_squares

@dataclass
class Spectrum:
    wavelength_nm: np.ndarray
    R: np.ndarray | None
    T: np.ndarray | None
    thickness_nm: float
    substrate_n: complex = 1.0
    ambient_n: complex = 1.0
    angle_deg: float = 0.0
    polarization: str = 's'
    sigma_R: float = .005
    sigma_T: float = .005
    name: str = 'spectrum'


def coherent_film(wavelength_nm, n_film, thickness_nm, substrate_n=1.,
                  ambient_n=1., angle_deg=0., polarization='s'):
    """Power R,T for coherent uniform film on semi-infinite substrate.
    For lossless ambient/substrate; absorbing external media not supported.
    """
    wl = np.asarray(wavelength_nm, dtype=float)
    nf = np.asarray(n_film, dtype=complex)
    n0, ns = complex(ambient_n), complex(substrate_n)
    if abs(n0.imag)>1e-12 or abs(ns.imag)>1e-12:
        raise ValueError('Ambient and semi-infinite substrate must be nonabsorbing')
    if np.any(np.real(nf)<=0) or np.any(np.imag(nf)<0):
        raise ValueError('Require n>0, k>=0')
    u=n0*np.sin(np.deg2rad(angle_deg))
    q0=np.sqrt(n0*n0-u*u+0j)
    qf=np.sqrt(nf*nf-u*u+0j)
    qs=np.sqrt(ns*ns-u*u+0j)
    # physical longitudinal-wave branch: positive real and positive imaginary
    qf=np.where(np.imag(qf)<0,-qf,qf)
    if polarization.lower() == 's':
        y0,yf,ys=q0,qf,qs
    elif polarization.lower() == 'p':
        y0,yf,ys=n0*n0/q0,nf*nf/qf,ns*ns/qs
    else:
        raise ValueError('polarization must be s or p')
    r01=(y0-yf)/(y0+yf)
    r12=(yf-ys)/(yf+ys)
    t01=2*y0/(y0+yf)
    t12=2*yf/(yf+ys)
    p=np.exp(2j*np.pi*qf*thickness_nm/wl)
    den=1+r01*r12*p*p
    r=(r01+r12*p*p)/den
    t=t01*t12*p/den
    R=np.abs(r)**2
    T=np.real(ys/y0)*np.abs(t)**2
    return R,T


def spline_design(wavelength_nm, lo, hi, n_basis=18, degree=3):
    x=np.asarray(wavelength_nm, dtype=float)
    if np.any((x<lo)|(x>hi)): raise ValueError('wavelength outside spline domain')
    internal=np.linspace(lo,hi,n_basis-degree+1)[1:-1]
    knots=np.r_[np.repeat(lo,degree+1),internal,np.repeat(hi,degree+1)]
    return BSpline(knots,np.eye(n_basis),degree)(x)


def fit_joint(spectra, n_basis=18, smoothness=1e-2, n_init=1.7, k_init=.02,
              max_nfev=300):
    """Joint fit, shared n(lambda), k(lambda), fixed sample thicknesses.
    Regularization is discrete second-difference, not a KK constraint.
    Returns result and callable optical_constants(wavelength_nm).
    """
    if not spectra: raise ValueError('At least one spectrum required')
    lo=min(np.min(s.wavelength_nm) for s in spectra)
    hi=max(np.max(s.wavelength_nm) for s in spectra)
    if lo>=hi: raise ValueError('Need multiple wavelengths')
    mats=[spline_design(s.wavelength_nm,lo,hi,n_basis) for s in spectra]
    def fun(x):
        an=x[:n_basis]; ak=x[n_basis:]
        out=[]
        for s,B in zip(spectra,mats):
            nc=B@an+1j*(B@ak)
            R,T=coherent_film(s.wavelength_nm,nc,s.thickness_nm,s.substrate_n,
                              s.ambient_n,s.angle_deg,s.polarization)
            if s.R is not None: out.extend((R-np.asarray(s.R))/s.sigma_R)
            if s.T is not None: out.extend((T-np.asarray(s.T))/s.sigma_T)
        out.extend(np.sqrt(smoothness)*np.diff(an,2))
        out.extend(np.sqrt(smoothness)*np.diff(ak,2))
        return np.asarray(out)
    x0=np.r_[np.full(n_basis,n_init),np.full(n_basis,k_init)]
    low=np.r_[np.full(n_basis,.5),np.zeros(n_basis)]
    high=np.r_[np.full(n_basis,6.),np.full(n_basis,5.)]
    result=least_squares(fun,x0,bounds=(low,high),max_nfev=max_nfev,
                         x_scale='jac',ftol=1e-10,xtol=1e-10,gtol=1e-10)
    def optical_constants(wavelength_nm):
        B=spline_design(wavelength_nm,lo,hi,n_basis)
        return B@result.x[:n_basis],B@result.x[n_basis:]
    return result,optical_constants


def local_sensitivity(wavelength_nm,n,k,thickness_nm,substrate_n=1.,
                      angle_deg=0.,polarization='s',step=1e-5):
    def f(nn,kk):
        r,t=coherent_film(np.array([wavelength_nm]),np.array([nn+1j*kk]),
                          thickness_nm,substrate_n,angle_deg=angle_deg,
                          polarization=polarization)
        return np.array([r[0],t[0]])
    dn=(f(n+step,k)-f(n-step,k))/(2*step)
    # central difference near k=0 is outside passive region, so one-sided there
    dk=((f(n,k+step)-f(n,k-step))/(2*step) if k>=step
        else (f(n,k+step)-f(n,k))/step)
    return np.column_stack([dn,dk])


def joint_information(configurations,wavelength_nm,n,k,sigma=.005):
    """Return singular values and covariance estimate for local n,k at wavelength.
    configurations is iterable of dicts containing thickness_nm and optional
    substrate_n, angle_deg, polarization. Assumes independent equal R/T noise.
    """
    J=np.vstack([local_sensitivity(wavelength_nm,n,k,**cfg) for cfg in configurations])
    U,s,Vh=np.linalg.svd(J/sigma,full_matrices=False)
    cov=np.linalg.pinv(J.T@J/sigma**2)
    return s,cov
