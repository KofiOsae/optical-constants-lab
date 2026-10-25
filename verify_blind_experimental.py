import numpy as np
import pandas as pd
from engine import FitSettings,fit_absorption_aware_blind
from optical_inversion_package.inversion import Spectrum

def test_blind_on_experimental_csv(path):
    raw=pd.read_csv(path)
    w=raw['wavelength_um'].to_numpy()*1000.
    l=w/1000.
    ns=np.sqrt(1+.6961663*l*l/(l*l-.0684043**2)+.4079426*l*l/(l*l-.1162414**2)+.8974794*l*l/(l*l-9.896161**2))
    s=Spectrum(w,raw['R'].to_numpy(),raw['T'].to_numpy(),2000.,substrate_n=ns)
    result=fit_absorption_aware_blind(s,FitSettings(n_basis=18,smoothness=.01,max_nfev=140),transparent_min_nm=680.)
    print({k:v for k,v in result.items() if not isinstance(v,np.ndarray)})
    return result

if __name__=='__main__':
    import sys
    test_blind_on_experimental_csv(sys.argv[1])
