import ast
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]

def test_reference_csv_wavelength_units():
    df=pd.read_csv(ROOT/'data'/'R6G_Singh.csv')
    assert list(df.columns)==['wavelength_um','n','k']
    assert len(df)>100
    assert np.isfinite(df[['wavelength_um','n','k']].to_numpy()).all()

def test_server_reference_handles_um_explicitly():
    source=(ROOT/'server.py').read_text()
    tree=ast.parse(source)
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='reference')
    code=compile(ast.Module(body=[fn],type_ignores=[]),str(ROOT/'server.py'),'exec')
    ns={'pd':pd,'np':np,'DATA':ROOT/'data'/'R6G_Singh.csv'}
    exec(code,ns)
    w,n,k=ns['reference']()
    assert np.isclose(w[0],193)
    assert np.isclose(w[-1],1690)
    assert len(w)==len(n)==len(k)
    assert np.all(np.diff(w)>0)
