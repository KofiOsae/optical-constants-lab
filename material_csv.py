"""Validated wavelength-dependent optical constants from literature or simulation CSV."""
import csv
import io
import numpy as np

def parse_nk_csv(content, wavelength_unit=None):
    """Return normalized arrays (wavelength_nm,n,k). Never guess ambiguous units."""
    if isinstance(content,bytes):
        content=content.decode('utf-8-sig')
    rows=[line for line in content.splitlines() if line.strip() and not line.lstrip().startswith('#')]
    if not rows:raise ValueError('Empty material file')
    dialect=csv.Sniffer().sniff('\n'.join(rows[:min(5,len(rows))]),delimiters=',;\t')
    reader=csv.DictReader(io.StringIO('\n'.join(rows)),dialect=dialect)
    if not reader.fieldnames:raise ValueError('Missing header')
    names={v.strip().lower().replace(' ','_'):v for v in reader.fieldnames}
    wkey=next((names[x] for x in ('wavelength_um','wavelength_nm','wl_um','wl_nm','lambda_um','lambda_nm') if x in names),None)
    if wkey is None:raise ValueError('Header must explicitly specify wavelength_um or wavelength_nm')
    unit='um' if wkey.strip().lower().endswith('_um') else 'nm'
    if wavelength_unit is not None and wavelength_unit!=unit:raise ValueError('Unit conflicts with header')
    if 'n' not in names or 'k' not in names:raise ValueError('Required columns n and k')
    data=[]
    for i,row in enumerate(reader,2):
        try: data.append([float(row[wkey]),float(row[names['n']]),float(row[names['k']])])
        except (ValueError,TypeError) as e:raise ValueError(f'Invalid numeric value at line {i}') from e
    a=np.asarray(data,dtype=float)
    if a.ndim!=2 or a.shape[0]<2:raise ValueError('At least two rows required')
    a[:,0]*=1000 if unit=='um' else 1
    if not np.all(np.isfinite(a)) or np.any(a[:,0]<=0) or np.any(a[:,1]<=0) or np.any(a[:,2]<0):
        raise ValueError('Require finite wavelength>0, n>0, k>=0')
    if np.any(np.diff(a[:,0])<=0):raise ValueError('Wavelengths must be strictly increasing and unique')
    return a[:,0],a[:,1],a[:,2]

def interpolate_nk(wavelength_nm,source_wavelength_nm,n,k):
    """Piecewise linear interpolation, strictly within reference wavelength domain."""
    x=np.asarray(wavelength_nm,float);w=np.asarray(source_wavelength_nm,float)
    if np.any(~np.isfinite(x)) or np.any(x<w[0]) or np.any(x>w[-1]):
        raise ValueError('Requested wavelength outside material reference range; extrapolation disabled')
    return np.interp(x,w,n)+1j*np.interp(x,w,k)
