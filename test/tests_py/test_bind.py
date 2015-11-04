import unittest
import random
import time
import numpy as np
from bspline_bind import *

class TestBind(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_non0_index(self):
        self.assertEqual((0, 6), non0_quad_index(0, 1, 3, 30))
        self.assertEqual((0, 6), non0_quad_index(0, 0, 3, 30))
        self.assertEqual((9, 12), non0_quad_index(2, 4, 3, 30))


    """
    def test_pdot(self):
        xs = np.array([random.random() for i in range(100000)])
        ys = np.array([random.random() for i in range(100000)])
        zs = np.array([random.random() for i in range(100000)])
        tm = time.clock()
        res = np.dot(xs, zs*ys)
        t0 = time.clock()
        res0 = np.dot(xs[100:1000], zs[100:1000]*ys[100:1000])
        t1 = time.clock()
        res1 = pdot3(xs, ys, zs, 100, 1000)
        t2 = time.clock()
        self.assertAlmostEqual(res0, res1)
        print "(dot, dot with slice, pdot)"
        print (t0-tm, t1-t0, t2-t1)
    """

    def test_outer_sum(self):
        b = np.array([1, 3, 4])
        a = np.array([2, -2])
        ans = np.array([3, 5, 6,
                        -1, 1, 2])
        calc = outer_sum(a, b)
        self.assertEqual(ans[0], calc[0])
        self.assertEqual(ans[1], calc[1])
        self.assertEqual(ans[2], calc[2])
        self.assertEqual(ans[3], calc[3])
        self.assertEqual(ans[4], calc[4])
        self.assertEqual(ans[5], calc[5])

        b = np.array(b, dtype=np.int32)
        a = np.array(a, dtype=np.int32)
        ans = np.array(ans, dtype=np.int32)
        calc = outer_sum(a, b)
        self.assertEqual(ans[0], calc[0])
        self.assertEqual(ans[1], calc[1])
        self.assertEqual(ans[2], calc[2])
        self.assertEqual(ans[3], calc[3])
        self.assertEqual(ans[4], calc[4])

if __name__ == '__main__':
    unittest.main()
