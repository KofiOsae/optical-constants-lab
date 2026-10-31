import io,json
import numpy as np
from server import app
def test_joint_fit_rejects_insufficient_substrate_ri_coverage():
    w=np.linspace(500,800,16)
    data=('wavelength_nm,R,T\n'+'\n'.join(f'{x},0.1,0.8' for x in w)).encode()
    cfg={'thickness_nm':180,'substrate_ri':{'wavelength_nm':[550,750],'n':[1.46,1.47],'k':[0,0]}}
    with app.test_client() as client:
        reply=client.post('/api/joint-fit',data={'spectra':(io.BytesIO(data),'test.csv'),'configurations':json.dumps([cfg]),'joint_method':'lorentz'},content_type='multipart/form-data')
    assert reply.status_code==400
    assert 'cover' in reply.get_json()['error'].lower()
def test_substrate_ri_interpolates_without_extrapolation():
    from multilayer_api import substrate_index
    cfg={'substrate_ri':{'wavelength_nm':[400,600,800],'n':[1.5,1.6,1.7],'k':[0,.01,.02]}}
    val=substrate_index(cfg,np.array([500.,700.]))
    assert np.allclose(val,[1.55+.005j,1.65+.015j])
