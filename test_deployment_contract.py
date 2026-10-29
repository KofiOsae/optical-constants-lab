from pathlib import Path
import ast
ROOT=Path(__file__).resolve().parent
def test_server_syntax():
    ast.parse((ROOT/'server.py').read_text())
def test_render_and_health_contract():
    assert 'gunicorn server:app' in (ROOT/'render.yaml').read_text()
    assert "@app.get('/health')" in (ROOT/'server.py').read_text()
def test_mobile_viewport():
    html=(ROOT/'web/index.html').read_text()
    assert 'name="viewport"' in html and '@media(max-width:800px)' in html

def test_browser_javascript_syntax():
    """Existing numerical tests cannot detect syntax errors in inline browser scripts."""
    import re, subprocess, tempfile
    for name in ('workspace.html','multilayer.html'):
        html=(ROOT/'web'/name).read_text(encoding='utf8')
        inline=re.findall(r'<script(?:\s[^>]*)?>(.*?)</script>',html,flags=re.S|re.I)
        for script in inline:
            if not script.strip():
                continue
            with tempfile.NamedTemporaryFile(mode='w',suffix='.js',encoding='utf8') as tmp:
                tmp.write(script)
                tmp.flush()
                check=subprocess.run(['node','--check',tmp.name],capture_output=True,text=True)
                assert check.returncode==0, f'{name}: {check.stderr}'

def test_multilayer_preview_um_csv_and_fit_range():
    import io, json
    import numpy as np
    from server import app
    client=app.test_client()
    rows=['wavelength_um,R,T']+[f'{w:.6f},0.2,0.7' for w in np.linspace(.4,1.,121)]
    payload='\n'.join(rows).encode()
    response=client.post('/api/multilayer/preview',data={
        'spectrum':(io.BytesIO(payload),'sample.csv')},content_type='multipart/form-data')
    assert response.status_code==200,response.get_json()
    obj=response.get_json()
    assert obj['point_count']==121
    assert np.isclose(obj['wavelength_nm'][0],400)
    assert np.isclose(obj['wavelength_nm'][-1],1000)
    assert np.allclose(obj['R'],.2)
    assert np.allclose(obj['T'],.7)
    response=client.post('/api/multilayer/joint-fit',data={
        'spectra':(io.BytesIO(payload),'sample.csv'),
        'configurations':json.dumps([{'thickness_nm':80,'substrate_n':1.46,
                                      'fit_range_nm':[405,425]}])},
        content_type='multipart/form-data')
    assert response.status_code==400,response.get_json()
    assert '8' in response.get_json()['error']


def test_multilayer_preview_requires_declared_wavelength_unit_and_two_channels():
    import io
    from server import app
    client=app.test_client()
    ambiguous='wavelength,R,T\n'+''.join(f'{w},0.2,0.7\n' for w in range(400,480,10))
    response=client.post('/api/multilayer/preview',data={
        'spectrum':(io.BytesIO(ambiguous.encode()),'sample.csv')},
        content_type='multipart/form-data')
    assert response.status_code==400
    incomplete='wavelength_nm,R\n'+''.join(f'{w},0.2\n' for w in range(400,480,10))
    response=client.post('/api/multilayer/preview',data={
        'spectrum':(io.BytesIO(incomplete.encode()),'sample.csv')},
        content_type='multipart/form-data')
    assert response.status_code==400
