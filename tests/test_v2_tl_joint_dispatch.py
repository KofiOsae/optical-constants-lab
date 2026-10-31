import io,json
import numpy as np
from server import app
from organic_dispersion import organic_nk
from optical_inversion_package.inversion import coherent_film

def test_tauc_lorentz_joint_dispatch():
    w=np.arange(450.,851.,25.)
    n,k=organic_nk(w,[2.35,18.,1.5,1.2,.32],'tauc_lorentz',1)
    R,T=coherent_film(w,n+1j*k,180.,1.46)
    csv=('wavelength_nm,R,T\n'+'\n'.join(f'{x},{r},{t}' for x,r,t in zip(w,R,T))).encode()
    config=[{'name':'tl','thickness_nm':180.,'substrate_n':1.46,'substrate_k':0}]
    with app.test_client() as client:
        response=client.post('/api/joint-fit',data={'spectra':(io.BytesIO(csv),'tl.csv'),'configurations':json.dumps(config),'joint_method':'tauc_lorentz','oscillators':'1','max_nfev':'6'},content_type='multipart/form-data')
    assert response.status_code==200,response.get_json()
    assert 'error' not in response.get_json()
