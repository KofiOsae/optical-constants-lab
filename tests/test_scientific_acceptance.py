"""Scientific acceptance separates forward-model parity from inverse accuracy."""
import numpy as np
from server import app
from historical_reference import archived_phase2

def test_archived_forward_model_parity_and_metrics():
    response=app.test_client().get('/api/benchmark/scientific-acceptance')
    assert response.status_code==200,response.get_json()
    data=response.get_json()
    assert data['passed'],data
    metrics=data['metrics']
    for key in ('R_train_rmse','T_train_rmse','R_holdout_rmse','T_holdout_rmse',
                'n_reference_mae','k_reference_mae','n_reference_rmse','k_reference_rmse'):
        assert np.isfinite(metrics[key]) and metrics[key]>=0
    assert metrics['forward_R_max_abs']<=data['limits']['forward_R_max_abs']
    assert metrics['forward_T_max_abs']<=data['limits']['forward_T_max_abs']
    assert data['reference_only'] is True
