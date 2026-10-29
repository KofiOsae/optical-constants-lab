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

def test_multilayer_api_distinguishes_fixed_and_fitted_thickness():
    import io, json
    import numpy as np
    from server import app
    from multilayer_forward import coherent_stack
    w=np.linspace(500.,750.,21)
    reflect,transmit=coherent_stack(w,[(1.8+.05j,100.)],1.5)
    rows=['wavelength_nm,R,T']+[f'{x:.8f},{r:.12f},{t:.12f}' for x,r,t in zip(w,reflect,transmit)]
    raw='\n'.join(rows).encode()
    cfg={'thickness_nm':100.,'thickness_bounds_nm':[90.,110.],
         'substrate_n':1.5}
    for enabled in (False,True):
        response=app.test_client().post('/api/multilayer/joint-fit',data={
            'spectra':(io.BytesIO(raw),'fixed.csv'),
            'configurations':json.dumps([cfg]),
            'model':'lorentz','oscillators':'1','max_nfev':'10',
            'fit_thickness':'true' if enabled else 'false'},
            content_type='multipart/form-data')
        assert response.status_code==200,response.get_json()
        body=response.get_json()
        p=body['predictions'][0]
        assert body['fit_thickness'] is enabled
        assert p['thickness_fitted'] is enabled
        assert p['nominal_thickness_nm']==100.
        if enabled:
            assert p['thickness_bounds_nm']==[90.,110.]
            assert 90.<=p['thickness_nm']<=110.
        else:
            assert p['thickness_nm']==100.
            assert p['thickness_bounds_nm'] is None

def test_multilayer_api_warm_start_reuse_and_domain_guard():
    import io, json
    import numpy as np
    from server import app
    from multilayer_forward import coherent_stack
    w=np.linspace(510.,800.,21)
    R,T=coherent_stack(w,[(1.85+.04j,120.)],1.5)
    payload=('wavelength_nm,R,T\n'+'\n'.join(f'{x:.6f},{r:.12f},{t:.12f}' for x,r,t in zip(w,R,T))).encode()
    config=json.dumps([{'thickness_nm':120.,'substrate_n':1.5}])
    client=app.test_client()
    def run(warm=None):
        data={'spectra':(io.BytesIO(payload),'sample.csv'),
              'configurations':config,'model':'spline',
              'n_basis':'4','max_nfev':'35'}
        if warm is not None:data['warm_start']=json.dumps(warm)
        return client.post('/api/multilayer/joint-fit',data=data,
                           content_type='multipart/form-data')
    first=run()
    assert first.status_code==200,first.get_json()
    saved=first.get_json()['warm_start_package']
    second=run(saved)
    assert second.status_code==200,second.get_json()
    assert second.get_json()['warm_start_used'] is True
    invalid={**saved,'wavelength_range_nm':[500.,800.]}
    reject=run(invalid)
    assert reject.status_code==400
    assert 'domain' in reject.get_json()['error'].lower()

def test_single_absorbing_film_joint_accepts_lorentz_and_direct():
    import io,json
    import numpy as np
    from server import app
    from optical_inversion_package.inversion import coherent_film
    w=np.linspace(500.,850.,21)
    n=1.8+0*w
    k=.08+0*w
    R,T=coherent_film(w,n+1j*k,135.,1.5)
    payload=('wavelength_nm,R,T\n'+'\n'.join(f'{x:.6f},{r:.12f},{t:.12f}' for x,r,t in zip(w,R,T))).encode()
    config=json.dumps([{'name':'absorbing','thickness_nm':135.,'substrate_n':1.5}])
    for method in ('direct','lorentz'):
        response=app.test_client().post('/api/joint-fit',data={
            'spectra':(io.BytesIO(payload),'absorbing.csv'),
            'configurations':config,'joint_method':method,
            'basis':'4','oscillators':'1','max_nfev':'25',
            'optimizer':'least_squares','causal_starts':'1'},
            content_type='multipart/form-data')
        assert response.status_code==200,(method,response.get_json())
        obj=response.get_json()
        assert len(obj['n'])>0 and len(obj['k'])>0
        assert len(obj['predictions'])==1


def test_joint_preview_and_tutorial_equations_are_browser_ready():
    import re
    page=(ROOT/'web'/'workspace.html').read_text()
    assert "async function previewJointUploads()" in page
    assert "/api/import-preview" in page
    assert "if(fs.length===1&&!['cauchy','sellmeier']" not in page
    equations=re.findall(r'<div class="equation">(.*?)</div>',page,re.S)
    assert len(equations)>=4
    assert all(('\\\\' not in expression) and ('\\[' in expression) and ('\\]' in expression) for expression in equations)


def test_multilayer_all_model_warm_start_frontend_controls():
    page=(ROOT/'web'/'multilayer.html').read_text()
    assert 'const warmCache=new Map()' in page
    assert "for(const id of ['model','basis','osc'])" in page
    assert "warmCache.get(warmKeyFromControls())" in page

def test_noise_copy_seed_reproducibility_and_bounds():
    import numpy as np
    import pytest
    from optical_inversion_package.inversion import Spectrum
    from custom_spectra import apply_rt_noise
    w=np.linspace(500.,700.,25)
    source=Spectrum(w,np.full(25,.2),np.full(25,.6),100.,substrate_n=1.5)
    first=apply_rt_noise(source,.002,1234,0)
    repeated=apply_rt_noise(source,.002,1234,0)
    independent=apply_rt_noise(source,.002,1234,1)
    unchanged=apply_rt_noise(source,0.,1234,0)
    assert np.array_equal(first.R,repeated.R)
    assert np.array_equal(first.T,repeated.T)
    assert not np.array_equal(first.R,independent.R)
    assert np.array_equal(unchanged.R,source.R)
    assert np.array_equal(unchanged.T,source.T)
    assert np.array_equal(source.R,np.full(25,.2))
    assert np.array_equal(source.T,np.full(25,.6))
    assert np.all((first.R>=0)&(first.R<=1))
    assert np.all((first.T>=0)&(first.T<=1))
    with pytest.raises(ValueError,match='Noise'):
        apply_rt_noise(source,.2,1234,0)


def test_joint_and_multilayer_noise_api_reproducible():
    import io,json
    import numpy as np
    from server import app
    from multilayer_forward import coherent_stack
    w=np.linspace(500.,800.,21)
    R,T=coherent_stack(w,[(1.9+.04j,120.)],1.5)
    csv=('wavelength_nm,R,T\n'+'\n'.join(f'{x:.6f},{r:.12f},{t:.12f}' for x,r,t in zip(w,R,T))).encode()
    client=app.test_client()
    def send(route):
        cfg=[{'thickness_nm':120.,'substrate_n':1.5}]
        d={'spectra':(io.BytesIO(csv),'sample.csv'),
           'configurations':json.dumps(cfg),
           'noise_sigma':'.002','noise_seed':'1234',
           'max_nfev':'12'}
        if route.endswith('joint-fit') and '/multilayer/' not in route:
            d.update({'joint_method':'lorentz','oscillators':'1','causal_starts':'1'})
        else:
            d.update({'model':'lorentz','oscillators':'1'})
        return client.post(route,data=d,content_type='multipart/form-data')
    for route in ('/api/joint-fit','/api/multilayer/joint-fit'):
        a,b=send(route),send(route)
        assert a.status_code==200,a.get_json()
        assert b.status_code==200,b.get_json()
        aa,bb=a.get_json(),b.get_json()
        assert aa['noise']['sigma']==.002 and aa['noise']['seed']==1234
        assert aa['predictions'][0]['R_measured']==bb['predictions'][0]['R_measured']
        assert aa['predictions'][0]['T_measured']==bb['predictions'][0]['T_measured']


def test_all_workflows_expose_noise_controls():
    root=ROOT/'web'
    joint=(root/'workspace.html').read_text()
    multi=(root/'multilayer.html').read_text()
    for html in (joint,multi):
        assert 'id="noise_sigma"' in html
        assert 'id="noise_seed"' in html
        assert "f.append('noise_sigma'" in html

def test_multilayer_staged_lorentz_api_high_to_full_wavelengths():
    import io,json
    import numpy as np
    from server import app
    from multilayer_forward import coherent_stack
    w=np.linspace(450.,950.,26)
    eps=2.5+.3*(2.1**2)/(2.1**2-(1239.841984/w)**2-1j*.2*(1239.841984/w))
    nk=np.sqrt(eps+0j)
    rr,tt=coherent_stack(w,[(1.65,60.),(nk,130.)],1.5)
    csv=('wavelength_nm,R,T\n'+'\n'.join(f'{x:.6f},{r:.12f},{t:.12f}' for x,r,t in zip(w,rr,tt))).encode()
    def post(lo=None,hi=None,previous=None):
        cfg={'thickness_nm':130.,'substrate_n':1.5,
             'before':[{'n':1.65,'thickness_nm':60.}]}
        if lo is not None:cfg['fit_range_nm']=[lo,hi]
        body={'spectra':(io.BytesIO(csv),'spectrum.csv'),
              'configurations':json.dumps([cfg]),'model':'lorentz',
              'oscillators':'1','max_nfev':'45'}
        if previous is not None:body['warm_start_spectrum']=json.dumps(previous)
        return app.test_client().post('/api/multilayer/joint-fit',data=body,
                                      content_type='multipart/form-data')
    first=post(690.,950.)
    assert first.status_code==200,first.get_json()
    stage=first.get_json()
    assert min(stage['wavelength_nm'])>=690
    assert not stage['warm_start_used']
    previous={'source':'previous_multilayer_fit','model':'lorentz',
              'wavelength_nm':stage['wavelength_nm'],
              'n':stage['n'],'k':stage['k']}
    second=post(previous=previous)
    assert second.status_code==200,second.get_json()
    extended=second.get_json()
    assert min(extended['wavelength_nm'])<min(stage['wavelength_nm'])
    assert extended['warm_start_used'] is True
    assert extended['warm_start_transferred'] is True
    assert extended['multistart']['attempts']>=1
    bad={**previous,'model':'spline'}
    rejected=post(previous=bad)
    assert rejected.status_code==400


def test_multilayer_preview_is_button_only_and_staged_controls_present():
    import re
    html=(ROOT/'web'/'multilayer.html').read_text(encoding='utf8')
    assert 'id="previewPanel" hidden' in html
    assert "$('previewButton').onclick=previewSpectra" in html
    assert "addEventListener('change',previewSpectra)" not in html
    assert "onclick=()=>runMultilayer('stage1')" in html
    assert "onclick=()=>runMultilayer('stage2')" in html
    assert "f.append('warm_start_spectrum'" in html
