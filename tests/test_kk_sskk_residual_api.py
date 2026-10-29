import numpy as np

def test_kk_sskk_dielectric_residual_export_from_multilayer_api():
    import io,json
    from server import app
    from multilayer_forward import coherent_stack
    w=np.linspace(500,900,31)
    R,T=coherent_stack(w,[(1.8+.03j,120.)],1.5)
    csv=('wavelength_nm,R,T\n'+'\n'.join(f'{wl},{r},{t}' for wl,r,t in zip(w,R,T))).encode()
    client=app.test_client()
    data={'spectra':(io.BytesIO(csv),'sample.csv'),
          'configurations':json.dumps([{'thickness_nm':120.,'substrate_n':1.5}]),
          'model':'lorentz','oscillators':'1','max_nfev':'18',
          'dispersion_diagnostic':'sskk','anchor_wavelength_nm':'700',
          'anchor_epsilon1':'2.8'}
    response=client.post('/api/multilayer/joint-fit',data=data,content_type='multipart/form-data')
    assert response.status_code==200,response.get_json()
    report=response.get_json()['dispersion_diagnostic']
    assert len(report['epsilon1_residual'])==len(report['wavelength_nm'])
    assert np.isfinite(report['rmse_full'])
    assert np.isfinite(report['rmse_interior'])
    assert report['interior_wavelength_nm'][0]>min(report['wavelength_nm'])
    assert report['interior_wavelength_nm'][1]<max(report['wavelength_nm'])
