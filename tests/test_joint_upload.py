"""End-to-end multi-spectrum API tests with withheld synthetic optical constants."""
import io,json
import numpy as np
from server import app
from optical_inversion_package.inversion import coherent_film

def synthetic_csv(d,angle=0,pol='s'):
    w=np.linspace(500,1100,61)
    n=1.75+.07*(800/w)**2
    k=.008+.035*np.exp(-((w-550)/120)**2)
    R,T=coherent_film(w,n+1j*k,d,1.46,angle_deg=angle,polarization=pol)
    rows=['wavelength_nm,R,T']+[f'{a:.8f},{b:.12f},{c:.12f}' for a,b,c in zip(w,R,T)]
    return io.BytesIO(('\n'.join(rows)+'\n').encode())

def test_joint_upload_two_thicknesses():
    client=app.test_client()
    configs=[dict(thickness_nm=500,substrate_n=1.46,angle_deg=0,polarization='s'),
             dict(thickness_nm=1000,substrate_n=1.46,angle_deg=0,polarization='s')]
    response=client.post('/api/joint-fit',data={
        'spectra':[(synthetic_csv(500),'thin.csv'),(synthetic_csv(1000),'thick.csv')],
        'configurations':json.dumps(configs),'basis':'8','max_nfev':'90',
        'smoothness':'0.01','optimizer':'least_squares'},content_type='multipart/form-data')
    assert response.status_code==200,response.get_json()
    data=response.get_json()
    assert len(data['predictions'])==2
    assert len(data['n'])==len(data['wavelength_nm'])
    assert np.all(np.isfinite(data['n'])) and np.all(np.isfinite(data['k']))
    assert data['weighted_rmse']<20

def test_joint_upload_rejects_mismatched_configs():
    response=app.test_client().post('/api/joint-fit',data={
        'spectra':[(synthetic_csv(500),'one.csv'),(synthetic_csv(1000),'two.csv')],
        'configurations':'[]'},content_type='multipart/form-data')
    assert response.status_code==400
