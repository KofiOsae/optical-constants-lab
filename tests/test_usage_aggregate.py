from server import app

def test_usage_summary_without_persistent_backend(monkeypatch):
    monkeypatch.delenv('USAGE_DATABASE_URL',raising=False)
    client=app.test_client()
    response=client.get('/api/usage/summary?days=30')
    assert response.status_code==200
    data=response.get_json()
    assert data['enabled'] is False
    assert 'not connected' in data['reason']
    assert client.post('/api/usage/event',json={'event':'visit'}).status_code==503

def test_event_rejects_personal_payload_and_invalid_names(monkeypatch):
    monkeypatch.delenv('USAGE_DATABASE_URL',raising=False)
    client=app.test_client()
    assert client.post('/api/usage/event',json={'event':'visit','ip':'1.2.3.4'}).status_code==400
    assert client.post('/api/usage/event',json={'event':'name'}).status_code==400
    assert client.get('/api/usage/summary?days=3').status_code==400

def test_only_successful_fits_are_eligible_for_recording(monkeypatch):
    import usage_aggregate
    seen=[]
    monkeypatch.setattr(usage_aggregate,'record',lambda event:seen.append(event) or True)
    with app.test_request_context('/api/joint-fit',method='POST'):
        from flask import Response
        usage_aggregate.count_completed_requests(Response(status=500))
        assert seen==[]
        usage_aggregate.count_completed_requests(Response(status=200))
    assert seen==['fit_joint']
