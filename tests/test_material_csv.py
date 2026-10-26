import unittest
import numpy as np
from material_csv import parse_nk_csv,interpolate_nk
class MaterialCsvTest(unittest.TestCase):
    def test_um_and_interpolation(self):
        w,n,k=parse_nk_csv('wavelength_um,n,k\n0.5,1.5,0\n0.6,1.7,0.2\n')
        np.testing.assert_allclose(w,[500,600])
        self.assertAlmostEqual(interpolate_nk([550],w,n,k)[0],1.6+.1j)
    def test_nm(self):
        w,n,k=parse_nk_csv('wavelength_nm,n,k\n500,1.5,0\n600,1.7,0.2')
        np.testing.assert_allclose(w,[500,600])
    def test_invalid(self):
        for data in ('wavelength_um,n,k\n0.5,1.5,0\n0.5,1.7,0.2',
                     'wavelength_um,n,k\n0.5,1.5,-1\n0.6,1.7,0',
                     'wavelength,n,k\n500,1.5,0\n600,1.7,0'):
            with self.assertRaises(ValueError):parse_nk_csv(data)
    def test_no_extrapolation(self):
        with self.assertRaises(ValueError):interpolate_nk([499],[500,600],[1.5,1.6],[0,0])
if __name__=='__main__':unittest.main()
