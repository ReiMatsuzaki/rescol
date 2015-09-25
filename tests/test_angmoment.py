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

    def test_y1mat_Yqk(self):
        self.assertAlmostEqual(
            -wigner.gaunt(2, 1, 3, -1, 0, 1, 10.0),
            y1mat_Yqk((2, 1), (1, 0), (3, 1)))

    def test_coupledY(self):
        y1 = CoupledY((1, 1), 2, 0)
        self.assertEqual(2, y1.L)
        self.assertEqual(1, y1.L1)

    def test_y2mat_YYq(self):
        y0 = CoupledY((0, 0), 0, 0)
        self.assertAlmostEqual(1.0, y2mat_Pq_r12(y0, 0, y0))
        self.assertAlmostEqual(0.0, y2mat_YYq(y0, 1, y0))
                               

if __name__ == '__main__':
    unittest.main()
