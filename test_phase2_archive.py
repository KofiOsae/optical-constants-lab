from historical_reference import archived_phase2

def test_original_phase2_coefficients_reproduce_archived_metrics():
    result=archived_phase2()
    assert result['coefficient_count']==50
    assert abs(result['thickness_nm']-999.3798755854683)<1e-6
    assert abs(result['train_rmse']-0.003662353448643851)<2e-8
    assert abs(result['holdout_rmse']-0.0029839159124723912)<2e-8

def test_archived_benchmark_api():
    from server import app
    response=app.test_client().get('/api/benchmark/phase2-archive')
    assert response.status_code==200
    assert response.json['reference_only'] is True
