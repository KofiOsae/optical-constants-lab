"""Multilayer forward-model regression tests. Run: python -m unittest tests.test_multilayer"""
import unittest
import numpy as np
from multilayer import multilayer_rt
from optical_inversion_package.inversion import coherent_film
class MultilayerTest(unittest.TestCase):
    def test_single_layer_matches_reference(self):
        w=np.linspace(450,1100,81)
        for pol in ('s','p'):
            for angle in (0.,30.):
                N=1.8+0.04j
                a,b=multilayer_rt(w,[{'n':N,'thickness_nm':350}],1.46,1.,angle,pol)
                c,d=coherent_film(w,np.full_like(w,N,dtype=complex),350,1.46,1.,angle,pol)
                np.testing.assert_allclose(a,c,rtol=1e-10,atol=1e-11)
                np.testing.assert_allclose(b,d,rtol=1e-10,atol=1e-11)
    def test_lossless_energy(self):
        w=np.linspace(450,1100,81)
        layers=[{'n':2.2,'thickness_nm':100},{'n':1.45,'thickness_nm':170},{'n':2.2,'thickness_nm':100}]
        for pol in ('s','p'):
            R,T=multilayer_rt(w,layers,1.46,1.,25.,pol)
            np.testing.assert_allclose(R+T,1.,atol=1e-10)
    def test_absorbing_stack(self):
        w=np.linspace(450,1100,81)
        R,T=multilayer_rt(w,[{'n':2.+.2j,'thickness_nm':100},{'n':1.5,'thickness_nm':150}],1.46)
        self.assertTrue(np.all(R>=0) and np.all(T>=0) and np.all(R+T<=1+1e-9))
if __name__=='__main__':unittest.main()
