"""Reproducible classical TMM demonstrations, independent of web hosting.

Run: python examples/classical_multilayer_demo.py
Produces CSVs in examples/output/ with wavelength_nm,R,T,A.
"""
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from multilayer import multilayer_rt
from scientific_diagnostics import dielectric

OUT=Path(__file__).resolve().parent/'output'
OUT.mkdir(exist_ok=True)
w=np.linspace(450,1100,131)
cases={
    'bare_glass': [],
    'single_lossless_film': [{'n':1.8,'thickness_nm':200}],
    'three_layer_dielectric': [{'n':2.2,'thickness_nm':100},{'n':1.45,'thickness_nm':170},{'n':2.2,'thickness_nm':100}],
    'two_layer_absorbing': [{'n':2.0+0.2j,'thickness_nm':100},{'n':1.5,'thickness_nm':150}],
}
N=np.sqrt(dielectric(w,[2.5,0.6,2.0,0.2],'lorentz',1)+0j)
cases['lorentz_film']=[{'n':N,'thickness_nm':500}]
for name,layers in cases.items():
    R,T=multilayer_rt(w,layers,substrate_n=1.46,ambient_n=1.0,angle_deg=0,polarization='s')
    A=1-R-T
    np.savetxt(OUT/(name+'.csv'),np.column_stack([w,R,T,A]),delimiter=',',header='wavelength_nm,R,T,A',comments='')
    if np.any(A < -1e-9):raise RuntimeError(f'Negative absorption in {name}')
    print(f'{name}: R=[{R.min():.6f},{R.max():.6f}], T=[{T.min():.6f},{T.max():.6f}], A=[{A.min():.6f},{A.max():.6f}]')
print('Outputs:',OUT)
