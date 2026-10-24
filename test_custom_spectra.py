import io
import numpy as np
from custom_spectra import parse_spectrum
from server import app

def test_csv_um_percent():
    rows='wavelength_um,R,T\n'+''.join(f'{v:.3f},20,70\n' for v in np.linspace(.4,.9,12))
    s=parse_spectrum('test.csv',rows.encode(),500,1.5,signal_unit='percent')
    assert abs(s.wavelength_nm[0]-400)<1e-8
    assert np.allclose(s.R,.2)
    assert np.allclose(s.T,.7)

def test_json_nm():
    import json
    data={'wavelength_nm':list(np.linspace(400,800,12)),'reflectance':[.2]*12,'transmittance':[.7]*12}
    s=parse_spectrum('test.json',json.dumps(data).encode(),700,1.0)
    assert len(s.wavelength_nm)==12

def test_upload_api():
    rows='wavelength_nm,R,T\n'+''.join(f'{v:.3f},0.2,0.7\n' for v in np.linspace(450,850,12))
    client=app.test_client()
    response=client.post('/api/custom-fit',data={'spectrum':(io.BytesIO(rows.encode()),'sample.csv'),'thickness_nm':'600','basis':'4','max_nfev':'5'},content_type='multipart/form-data')
    assert response.status_code==200,response.json
    assert len(response.json['n'])==12
    assert 'n_true' not in response.json
