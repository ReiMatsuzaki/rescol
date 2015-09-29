import unittest
import time
import scipy.sparse.linalg as la
from bspline import *
from hmole import *
from angmoment import triangle_q, CoupledY


class TestHMole(unittest.TestCase):
    def setUp(self):
        pass

    def test_h_mole_plus_ground(self):
        order = 4
        num = 11
        rmax = 30.0
        r_h_mole = 1.4
        basis_set = BSplineSet(order, lin_knots(0.0, rmax, num))
        (h, s) = mat_h2_plus(r_h_mole, basis_set, [0, 2, 4, 6, 8])
        (es, vecs) = la.eigsh(h, 5, s, sigma=-2.0)
        self.assertTrue(abs(es[0] - (-1.1466859)) < 0.01)

    def test_h_mole_ground(self):
        order = 4
        num_knots = 11
        rmax = 30.0
        """
        order = 6
        num_knots = 31
        rmax = 30.0
        """
        bspline_set = BSplineSet(order, lin_knots(0.0, rmax, num_knots))
        y_list = [CoupledY((L1, L2), L, 0) for L in [0, 2, 4]
                  for L1 in [0, 2] for L2 in [0, 2]
                  if triangle_q(L1, L2, L)]
        bond_length = 1.4
        (h, s) = mat_h2(bond_length, bspline_set, y_list)
        (es, vecs) = la.eigsh(h, 5, s, sigma=-3.0)
        # 1.675 is come from 2015/9/28_rescol_check_h2mole.ipynb
        self.assertTrue(abs(es[0] - (-1.675)) < 0.01)


if __name__ == '__main__':
    unittest.main()
