import unittest
import scipy.sparse.linalg as la
from bspline import *
from hmole import *


class TestHMole(unittest.TestCase):
    def setUp(self):
        pass

    def test_h_mole_plus_ground(self):
        order = 8
        lmax = 8
        rmax = 60.0
        num = 140
        r_h_mole = 2.0
        basis_set = BSplineSet(order, lin_knots(0.0, rmax, num))
        (h, s) = mat_h2_plus(r_h_mole, basis_set, [0, 2, 4, 6, 8])
        (es, vecs) = la.eigsh(h, 5, s, sigma=-2.0)
        self.assertTrue(abs(es[0] - (-1.1026)) < 0.01)

if __name__ == '__main__':
    unittest.main()
