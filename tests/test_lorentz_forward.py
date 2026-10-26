"""Deterministic Lorentz forward-model smoke tests; run: python -m unittest tests.test_lorentz_forward"""
import unittest
import numpy as np
from scientific_diagnostics import dielectric
from engine import coherent_film

class LorentzForwardTest(unittest.TestCase):
    def test_passive_resonance(self):
        w=np.linspace(500,1100,61)
        p=[2.5,.6,2.0,.2]
        eps=dielectric(w,p,'lorentz',1)
        self.assertTrue(np.all(np.isfinite(eps)))
        self.assertTrue(np.all(np.imag(eps)>=0))
        N=np.sqrt(eps.astype(complex))
        self.assertTrue(np.all(np.imag(N)>=0))
        self.assertGreater(np.max(np.imag(N)),.01)
    def test_forward_energy(self):
        w=np.linspace(500,1100,61)
        N=np.sqrt(dielectric(w,[2.5,.6,2.,.2],'lorentz',1).astype(complex))
        for d in (500.,1000.,2000.):
            R,T=coherent_film(w,N,d,1.46,1.0,0.0,'s')
            self.assertTrue(np.all(np.isfinite(R)) and np.all(np.isfinite(T)))
            self.assertTrue(np.all(R>=-1e-8) and np.all(T>=-1e-8))
            self.assertTrue(np.all(R+T<=1+1e-7))
if __name__=='__main__':
    unittest.main()
