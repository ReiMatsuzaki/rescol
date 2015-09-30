import unittest
import scipy.sparse.linalg as la
from fdiff import *


class TestFiniteDiff(unittest.TestCase):

    def setUp(self):
        self.fdiff = FDiffSet(10, 20.0)

    def test_FDiffSet(self):

        self.assertEqual(10, len(self.fdiff.xs))
        self.assertEqual(2.0, self.fdiff.xs[0])
        self.assertEqual(20.0, self.fdiff.xs[-1])
        self.assertAlmostEqual(2.0, self.fdiff.h)

    def test_vmat(self):
        vmat = self.fdiff.v_mat(lambda x: 1.0/x)
        self.assertEqual((10, 10), vmat.shape)
        self.assertAlmostEqual(1.0/4.0, vmat.toarray()[1, 1])

    def test_d2mat(self):
        d2mat = self.fdiff.d2_mat()
        self.assertEqual((10, 10), d2mat.shape)
        self.assertAlmostEqual(-2/(2.0*2.0), d2mat.toarray()[0, 0])
        self.assertAlmostEqual(1/(2.0*2.0), d2mat.toarray()[0, 1])

    def test_hatom(self):
        basis_set = FDiffSet(1000, 30.0)
        d2_mat = basis_set.d2_mat()
        v_mat = -basis_set.en_mat(0, 0.0)
        h_mat = -0.5*d2_mat + v_mat
        (es, vecs) = la.eigsh(h_mat, 2, sigma=-0.6)
        self.assertTrue(abs(-0.5-es[0]) < 0.001, "es[0]={0}".format(es[0]))
        self.assertTrue(abs(-0.125-es[1]) < 0.001, "es[1]={0}".format(es[1]))


if __name__ == '__main__':
    unittest.main()
