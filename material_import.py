"""Shared, network-independent spectral table ingestion and inspection."""
import io
import json
import re
import numpy as np
import pandas as pd

MAX_BYTES = 1_000_000

def read_table(filename, content):
    if not content or len(content) > MAX_BYTES:
        raise ValueError('File must contain data and be no larger than 1 MB')
    filename = (filename or '').lower()
    if filename.endswith('.json'):
        obj = json.loads(content.decode('utf-8-sig'))
        if isinstance(obj, dict) and isinstance(obj.get('data'), list):
            obj = obj['data']
        if not isinstance(obj, (list, dict)):
            raise ValueError('JSON must contain records or arrays of columns')
        df = pd.DataFrame(obj)
    elif filename.endswith(('.csv', '.tsv', '.txt')):
        df = pd.read_csv(io.BytesIO(content), sep=None, engine='python', comment='#')
    else:
        raise ValueError('Supported formats: CSV, TSV, TXT, JSON')
    if not 2 <= len(df) <= 20000:
        raise ValueError('Require 2–20,000 data points')
    return df

def normalize(column):
    return re.sub('[^a-z0-9]', '', str(column).lower())

def inspect_material(filename, content, wavelength_unit='auto', signal_unit='auto', column_map=None):
    df = read_table(filename, content)
    normalized = {normalize(c): c for c in df.columns}
    aliases = {
        'wavelength': ['wavelengthnm', 'wavelengthum', 'wavelength', 'wl', 'lambda', 'lambdanm', 'lambdaum'],
        'n': ['n', 'refractiveindex', 'indexn'],
        'k': ['k', 'extinctioncoefficient', 'extinction', 'indexk'],
        'R': ['r', 'reflectance', 'reflection', 'rpercent', 'reflectancepercent', 'reflectionpercent'],
        'T': ['t', 'transmittance', 'transmission', 'tpercent', 'transmittancepercent', 'transmissionpercent'],
        'absorbance': ['absorbance', 'opticaldensity', 'od'],
        'absorptance': ['absorptance', 'absorptionfraction'],
        'psi': ['psi', 'psideg', 'ellipsometricpsi'],
        'delta': ['delta', 'deltadeg', 'ellipsometricdelta'],
        'eps1': ['eps1', 'epsilon1', 'realpermittivity'],
        'eps2': ['eps2', 'epsilon2', 'imagpermittivity'],
    }
    column_map = column_map or {}
    resolved = {}
    for name, keys in aliases.items():
        specified = column_map.get(name)
        if specified is not None:
            if specified not in df.columns:
                raise ValueError(f'Unknown column for {name}: {specified}')
            resolved[name] = specified
        else:
            resolved[name] = next((normalized[k] for k in keys if k in normalized), None)
    if not resolved['wavelength']:
        raise ValueError('Wavelength column not found; specify column mapping')
    w = pd.to_numeric(df[resolved['wavelength']], errors='coerce').to_numpy(dtype=float, copy=True)
    if not np.all(np.isfinite(w)):
        raise ValueError('Wavelength contains invalid values')
    if wavelength_unit == 'auto':
        key = normalize(resolved['wavelength'])
        if 'um' in key:
            wavelength_unit = 'um'
        elif 'nm' in key:
            wavelength_unit = 'nm'
        else:
            raise ValueError('Ambiguous wavelength unit: choose nm or µm explicitly')
    if wavelength_unit not in ('nm', 'um'):
        raise ValueError('Wavelength unit must be nm or um')
    if wavelength_unit == 'um':
        w *= 1000
    order = np.argsort(w)
    w = w[order]
    if w[0] <= 0 or np.any(np.diff(w) <= 0):
        raise ValueError('Wavelengths must be positive and unique')
    result = {'filename': filename, 'wavelength_nm': w.tolist(), 'columns': {}, 'mapping': {}, 'warnings': []}
    for key, column in resolved.items():
        if key == 'wavelength' or column is None:
            continue
        arr = pd.to_numeric(df[column], errors='coerce').to_numpy(dtype=float)[order]
        if not np.all(np.isfinite(arr)):
            raise ValueError(f'{column} contains missing or invalid values')
        if key in ('R', 'T', 'absorptance'):
            unit = signal_unit
            if unit == 'auto':
                unit = 'percent' if ('percent' in str(column).lower() or '%' in str(column)) else 'fraction'
            if unit not in ('fraction', 'percent'):
                raise ValueError('Signal unit must be fraction or percent')
            if unit == 'percent':
                arr /= 100
            if np.any((arr < 0) | (arr > 1)):
                raise ValueError(f'{column} must lie between 0 and 1 after conversion; check units')
        if key == 'k' and np.any(arr < 0):
            raise ValueError('Extinction coefficient k must be nonnegative for passive materials')
        result['columns'][key] = arr.tolist()
        result['mapping'][key] = str(column)
    cols = result['columns']
    if not cols:
        raise ValueError('No supported optical measurements or material constants found')
    if 'R' in cols and 'T' in cols and np.any(np.asarray(cols['R']) + np.asarray(cols['T']) > 1.02):
        raise ValueError('R+T exceeds 1.02; verify units and calibration')
    if 'n' in cols and 'k' in cols:
        result['kind'] = 'optical_constants'
    elif 'R' in cols or 'T' in cols:
        result['kind'] = 'film_spectrum'
        result['warnings'].append('A single film spectrum may not uniquely identify n, k and thickness.')
    elif 'psi' in cols and 'delta' in cols:
        result['kind'] = 'ellipsometry'
        result['warnings'].append('Ellipsometry inverse fitting is not implemented in this version; import and preview only.')
    else:
        result['kind'] = 'other_optical_data'
        result['warnings'].append('Imported for inspection; no compatible inverse solver is implemented for these columns.')
    result['point_count'] = len(w)
    result['wavelength_range_nm'] = [float(w[0]), float(w[-1])]
    return result
