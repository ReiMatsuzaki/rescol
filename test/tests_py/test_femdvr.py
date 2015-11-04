import unittest
import numpy as np
from scipy.sparse import coo_matrix, bmat
from femdvr import *


class TestUtils(unittest.TestCase):
    def setUp(self):
        pass

    def test_lin_knots(self):
        ts = lin_knots(1.0, 4.0, 4)  # => [1, 2, 3, 4]
        self.assertEqual(4, len(ts))
        self.assertAlmostEqual(1.0, ts[0])
        self.assertAlmostEqual(2.0, ts[1])

    def test_lobgauss(self):
        (xs, ws) = lobgauss(4)
        self.assertEqual(4, len(xs))
        self.assertEqual(4, len(ws))
        self.assertAlmostEqual(-1.0, xs[0])
        self.assertAlmostEqual(1.0/6.0, ws[0])
        self.assertAlmostEqual(-1.0/np.sqrt(5), xs[1])
        self.assertAlmostEqual(5.0/6.0, ws[1])
        self.assertAlmostEqual(+1.0/np.sqrt(5), xs[2])
        self.assertAlmostEqual(5.0/6.0, ws[2])
        self.assertAlmostEqual(+1.0, xs[3])
        self.assertAlmostEqual(1.0/6.0, ws[3])

    def test_lobgauss2(self):
        (xs, ws) = lobgauss(4, 0.0, 2.0)
        self.assertAlmostEqual(0.0, xs[0])
        self.assertAlmostEqual(1.0/6.0, ws[0])

    def test_f_to_chi(self):
        # see project/rescol/calc/fem_dvr.ipynb in elnote
        n = 4   # number of quad points in each elements
        ne = 3  # number of element
        A = coo_matrix(np.array([[i/n+ne*(j/n)
                                  for i in range(n*ne)]
                                 for j in range(n*ne)]))
        B = f_to_chi(A, n, ne)
        npB = B.toarray()
        self.assertEqual((8, 8), B.shape)
        self.assertAlmostEqual(0.0, npB[0, 0])
        self.assertAlmostEqual(5.0, npB[2, 3])
        self.assertAlmostEqual(20., npB[5, 2])


class _TestFemDvrSet(unittest.TestCase):
    def setUp(self):
        ts = lin_knots(0.0, 10.0, 9)
        n = 4
        self.n = n
        self.ne = len(ts) - 1
        self.femdvr_set = FemDvrSet(n, ts)

    def test_femdvr_s_fmat(self):
        n = self.n
        ne = self.ne
        s_mat = self.femdvr_set.s_fmat()
        self.assertEqual((n*ne, n*ne), s_mat.shape)
        self.assertEqual(n*ne, len(s_mat.data[0]),
                         "{0}, {1}\n {2}"
                         .format(n*ne, len(s_mat.data), s_mat.data))

    def test_femdvr_v_fmat(self):
        n = self.n
        ne = self.ne
        v_mat = self.femdvr_set.v_fmat(lambda r: 1.0/(r+0.00000001))
        self.assertEqual((n*ne, n*ne), v_mat.shape)
        self.assertEqual(n*ne, len(v_mat.data[0]))

    def test_femdvr_d2_fmat(self):
        n = self.n
        ne = self.ne
        d2_mat = self.femdvr_set.d2_fmat()
        self.assertEqual((n*ne, n*ne), d2_mat.shape)
        self.assertEqual(n**2*ne, len(d2_mat.data))


if __name__ == '__main__':
    unittest.main()
