import numpy as np
from engine import phase2_hybrid_seed

def test_hybrid_seed_uses_global_below_690_and_phase1_above():
    w=np.linspace(300,1000,201)
    global_n=np.full(len(w),1.65)
    global_k=np.full(len(w),.03)
    p=dict(a=1.7,b=0.,c=0.,q0=np.log(.002),q1=0.)
    q=phase2_hybrid_seed(w,w,global_n,global_k,p)
    kn=np.linspace(300,1000,25)
    low=kn<=690
    assert np.allclose(q[:25][low],1.65)
    assert np.allclose(np.exp(q[25:][low]),.03)
    assert np.allclose(q[:25][~low],1.7)
    assert np.allclose(np.exp(q[25:][~low]),.002)
