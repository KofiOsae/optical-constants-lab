import unittest
import numpy as np
from optical_inversion_package.inversion import coherent_film

class SubstrateDispersionTests(unittest.TestCase):
    def test_scalar_and_constant_array_equivalent(self):
        w=np.linspace(350.,1200.,51)
        z=1.8+.02j+np.zeros_like(w)
        a=coherent_film(w,z,2000.,substrate_n=1.46)
        b=coherent_film(w,z,2000.,substrate_n=np.full_like(w,1.46))
        for x,y in zip(a,b):
            np.testing.assert_allclose(x,y,rtol=1e-12,atol=1e-12)

    def test_lossless_energy_conservation_dispersive_silica(self):
        w=np.linspace(300.,1400.,201)
        l=w/1000.
        ns=np.sqrt(1+.6961663*l*l/(l*l-.0684043**2)+
                   .4079426*l*l/(l*l-.1162414**2)+
                   .8974794*l*l/(l*l-9.896161**2))
        r,t=coherent_film(w,np.full_like(w,1.8,dtype=complex),2000.,substrate_n=ns)
        np.testing.assert_allclose(r+t,1.,rtol=0,atol=1e-11)

    def test_reject_wrong_dispersion_shape(self):
        with self.assertRaises(ValueError):
            coherent_film(np.linspace(400.,800.,10),np.full(10,1.8+.01j),2000.,
                          substrate_n=np.full(9,1.46))

if __name__=='__main__':
    unittest.main()
