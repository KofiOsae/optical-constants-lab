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
