"""Number-weighted radius distributions for dilute sphere ensembles."""
import numpy as np
from scipy.special import ndtr,gammainc
from scipy.stats import weibull_min

def radius_bins(kind, characteristic_nm, spread=0.2, count=61,
                second_radius_nm=None, mixture_fraction=0.5,
                histogram_radii_nm=None,histogram_number_weights=None):
    if kind=="histogram":
        r=np.asarray(histogram_radii_nm,float);w=np.asarray(histogram_number_weights,float)
        if r.ndim!=1 or w.shape!=r.shape or not len(r) or np.any(~np.isfinite(r)) or np.any(r<=0) or np.any(~np.isfinite(w)) or np.any(w<0) or w.sum()<=0:
            raise ValueError("Histogram requires positive radii and nonnegative number weights")
        return r,w/w.sum()
    r0=float(characteristic_nm);s=float(spread)
    if not np.isfinite(r0) or r0<=0 or not np.isfinite(s) or s<0 or not 3<=count<=301:
        raise ValueError("Invalid radius, spread or bin count")
    if kind=="single" or s==0:return np.array([r0]),np.array([1.])
    if kind=="lognormal":
        from nanoparticle_optics.ensemble import lognormal_radius_bins
        return lognormal_radius_bins(r0,max(1.,1.+s),count=count)
    if kind=="bimodal":
        r2=float(second_radius_nm if second_radius_nm is not None else 2*r0)
        frac=float(mixture_fraction)
        if r2<=0 or not np.isfinite(r2) or not 0<=frac<=1:raise ValueError("Invalid bimodal parameters")
        a,wa=radius_bins("lognormal",r0,s,count)
        b,wb=radius_bins("lognormal",r2,s,count)
        return np.r_[a,b],np.r_[frac*wa,(1-frac)*wb]
    if kind not in ("normal","gamma","weibull"):raise ValueError("Unknown distribution")
    # Quantile-spaced bins avoid negative Gaussian radii and preserve probability mass.
    from scipy.stats import norm,gamma
    q=(np.arange(count)+.5)/count
    if kind=="normal":
        # Truncated normal on r>0, with characteristic_nm the parent mean.
        from scipy.stats import truncnorm
        sigma=r0*s
        r=truncnorm.ppf(q,-r0/sigma,np.inf,loc=r0,scale=sigma)
    elif kind=="gamma":
        shape=1/s**2;r=gamma.ppf(q,a=shape,scale=r0/shape)
    else:
        # spread is inverse Weibull shape; characteristic is scale.
        r=weibull_min.ppf(q,c=1/s,scale=r0)
    if not np.all(np.isfinite(r)) or np.any(r<=0):raise ValueError("Distribution numerical failure")
    return r,np.full(count,1/count)
