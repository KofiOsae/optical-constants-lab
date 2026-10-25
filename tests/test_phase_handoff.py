import io
import json
import numpy as np
from server import app
from optical_inversion_package.inversion import coherent_film

def test_phase1_thickness_range_tracks_specimen():
    from phase1_fresnel import phase1_fit
    w=np.linspace(850,1400,22)
    n=1.68+0.03*(800/w)**2
    k=np.full_like(w,0.001)
    R,T=coherent_film(w,n+1j*k,2500,1.46)
    out=phase1_fit(w,R,T,max_nfev=25,substrate_n=1.46,thickness_estimate_nm=2500)
    assert 1250<=out['thickness_nm']<=3750

def test_phase2_rejects_stale_phase1_approval():
    w=np.linspace(800,1400,12)
    R,T=coherent_film(w,np.full_like(w,1.7)+.001j,1000,1.46)
    data='wavelength_nm,R,T\n'+'\n'.join(f'{a},{b},{c}' for a,b,c in zip(w,R,T))
    payload={'spectrum':(io.BytesIO(data.encode()),'sample.csv'),'thickness_nm':'1000',
             'substrate_n':'1.46','reconstruction_model':'historical_logk',
             'wavelength_min_nm':'800','wavelength_max_nm':'1400',
             'approved_phase1':json.dumps({'parameters':{'a':1.7,'b':0,'c':0,'q0':-7,'q1':0,'d_nm':1000},
                                           'wavelength_nm':w.tolist(),'R_measured':[0]*12,'T_measured':[0]*12,
                                           'substrate_n':1.46,'structure':'supported'})}
    response=app.test_client().post('/api/custom-fit',data=payload,content_type='multipart/form-data')
    assert response.status_code==400
    assert 'Phase I input data changed' in response.get_json()['error']
