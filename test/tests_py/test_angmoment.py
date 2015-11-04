import unittest
import random
import time
import numpy as np
import sympy.physics.wigner as wigner
from angmoment import *


class TestAngMoment(unittest.TestCase):

    def setUp(self):
        pass

    def test_prod_LLL(self):
        self.assertAlmostEqual((2+1)*(2*3+1)*(2*4+1),
                               prod_LLL(1, 3, 4))

    def test_triangle_q(self):
        self.assertTrue(triangle_q(2, 2, 2))
        self.assertTrue(triangle_q(0, 2, 2))
        self.assertTrue(triangle_q(2, 0, 2))
        self.assertTrue(triangle_q(2, 2, 0))
        self.assertTrue(triangle_q(2, 1, 1))
        self.assertTrue(triangle_q(1, 2, 1))
        self.assertTrue(triangle_q(1, 1, 2))
        self.assertFalse(triangle_q(3, 1, 1))
        self.assertFalse(triangle_q(1, 0, 0))
        self.assertFalse(triangle_q(2, 0, 1))

    def test_zero_YYY_q(self):
        self.assertTrue(zero_YYY_q(1, 0, 0))
        self.assertTrue(zero_YYY_q(3, 3, 3))

    def test_non0_ls(self):
        self.assertEqual([3], ls_non_zero_YYY(0, 3))
        self.assertEqual([1, 3, 5], ls_non_zero_YYY(2, 3))
        self.assertEqual([0, 2, 4, 6], ls_non_zero_YYY(3, 3))

    def test_y1mat_Yqk(self):
        self.assertAlmostEqual(
            (-1)**(2+1)*wigner.gaunt(2, 1, 3, -1, 0, 1, 10.0),
            y1mat_Yqk((2, 1), (1, 0), (3, 1)))
        yyy = y1mat_Yqk((2, 1), (1, 0), (3, 1))
        self.assertEqual(type(yyy), np.float64)

    def test_y1mat_Pq(self):
        self.assertAlmostEqual(1.0, y1mat_Pq((0, 0), 0, (0, 0)))
        self.assertAlmostEqual(1.0, y1mat_Pq((1, 0), 0, (1, 0)))

    def test_coupledY(self):
        y1 = CoupledY((1, 1), 2, 0)
        self.assertEqual(2, y1.L)
        self.assertEqual(1, y1.L1)

    def test_coupledY_exchange(self):
        y1 = CoupledY((2, 1), 3, 0)
        y2 = y1.exchange()
        self.assertEqual(1, y2.L1)
        self.assertEqual(2, y1.L1)

    def test_get_coupledY(self):
        # y[(0,0), 0, 0], y[(1,1), 1, 0]
        # y[(2,0), 2, 0], y[(1,0), 1, 0], y[(0,2), 2, 0]
        ys = get_coupledY_set(2, 0, True, True)
        self.assertEqual(5, len(ys))

    def test_y2mat_YYq(self):
        y0 = CoupledY((0, 0), 0, 0)
        self.assertAlmostEqual(1.0, y2mat_Pq_r12(y0, 0, y0))
        self.assertAlmostEqual(0.0, y2mat_YYq(y0, 1, y0))

    def test_Y2mat_YYq_slow(self):
        yp = CoupledY((1, 1), 2, 0)
        y = CoupledY((3, 1), 2, 0)
        self.assertAlmostEqual(
            y2mat_YYq_slow(yp, 2, y),
            y2mat_YYq(yp, 2, y))

    def test_y2mat_Y1qk(self):
        yp = CoupledY((1, 1), 1, 0)
        y = CoupledY((3, 1), 3, 0)
        self.assertTrue(
            abs(y2mat_Y1qk(yp, 2, 0, y)) > 0.00001)
        self.assertAlmostEqual(
            y2mat_Y1qk_slow(yp, 2, 0, y),
            y2mat_Y1qk(yp, 2, 0, y))

if __name__ == '__main__':
    unittest.main()
