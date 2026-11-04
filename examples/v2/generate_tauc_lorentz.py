"""Generate an explicitly synthetic Tauc–Lorentz single-film example.

Usage: python examples/v2/generate_tauc_lorentz.py
Outputs RT CSV and known n,k reference CSV; the latter is reference model truth,
not independently measured values. Intended for v2 inverse-recovery experiments.
"""
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
import numpy as np
from organic_dispersion import organic_nk
from multilayer_forward import coherent_stack

ROOT=Path(__file__).resolve().parent
W=np.arange(380.,951.,5.)
PARAMS=[2.35,18.0,1.50,1.20,.32]  # [epsilon_inf,A,Eg,E0-Eg,C]
N,K=organic_nk(W,PARAMS,'tauc_lorentz',1)
R,T=coherent_stack(W,[(N+1j*K,180.)],substrate_n=1.46,ambient_n=1.)
np.savetxt(ROOT/'tauc_lorentz_180nm.csv',np.column_stack((W,R,T)),
           delimiter=',',header='wavelength_nm,R,T',comments='',fmt='%.10g')
np.savetxt(ROOT/'tauc_lorentz_180nm_reference_nk.csv',np.column_stack((W,N,K)),
           delimiter=',',header='wavelength_nm,n,k',comments='',fmt='%.10g')
print('Generated synthetic TL reference; max R+T:',float(np.max(R+T)))
