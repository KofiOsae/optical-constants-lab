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
        inline=re.findall(r'<script(?:\\s[^>]*)?>(.*?)</script>',html,flags=re.S|re.I)
        for script in inline:
            if not script.strip():
                continue
            with tempfile.NamedTemporaryFile(mode='w',suffix='.js',encoding='utf8') as tmp:
                tmp.write(script)
                tmp.flush()
                check=subprocess.run(['node','--check',tmp.name],capture_output=True,text=True)
                assert check.returncode==0, f'{name}: {check.stderr}'
