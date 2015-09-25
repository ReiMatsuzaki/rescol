import unittest
import time
import scipy.sparse.linalg as la
import scipy.sparse
from bspline_bind import *
from bspline import *


class TestBSpline(unittest.TestCase):

    def setUp(self):
        self.knots = make_knots((0.0, 3), 1.1, 2.2, 3.3, (4.5, 3))

    def test_make_knots(self):
        knots = make_knots(0.0, (3.0, 3), (2.0, 2))
        self.assertAlmostEqual(0.0, knots[0][0])
        self.assertEqual(1, knots[0][1])
        self.assertAlmostEqual(3.0, knots[1][0])
        self.assertEqual(3, knots[1][1])

    def test_lin_knots(self):
        knots = lin_knots(0.0, 5.0, 6)(3)
        ref = make_knots((0.0, 3), 1.0, 2.0, 3.0, 4.0, (5.0, 3))
        self.assertEqual(ref, knots)

    def test_intervals(self):
        intervals = make_intervals(self.knots)
        self.assertAlmostEqual(0.0, intervals[0][0])
        for i in range(3):
            self.assertAlmostEqual(intervals[0][1], intervals[1][0])
        self.assertAlmostEqual(4.5, intervals[3][1])
        self.assertAlmostEqual(3.3, intervals[3][0])

    def test_make_ts(self):
        knots = [(1.0, 3), (2.0, 1), (2.5, 2), (4.0, 3)]
        flattened = make_ts(knots)
        self.assertEqual(3+1+2+3, len(flattened))
        self.assertAlmostEqual(1.0, flattened[0])
        self.assertAlmostEqual(4.0, flattened[-1])

    def test_calc_bspline(self):
        # see paper
        knots = make_knots((0.0, 3), 1.0, 2.0, 3.0, 4.0, (5.0, 3))
        self.assertAlmostEqual(1.0, calc_bspline(3, knots, 0, 0.0))
        self.assertAlmostEqual(0.0, calc_bspline(3, knots, 0, 1.0))
        self.assertAlmostEqual(0.0, calc_bspline(3, knots, 6, 4.0))
        x = 0.34
        self.assertAlmostEqual(0.5*x*x, calc_bspline(3, knots, 2, x))
        self.assertAlmostEqual(x, calc_deriv_bspline(3, knots, 2, x))
        x = 2.44
        self.assertAlmostEqual(0.5*x*x-3*x+9.0/2.0,
                               calc_bspline(3, knots, 2, x))
        self.assertAlmostEqual(x-3, calc_deriv_bspline(3, knots, 2, x))
        x = 3.44
        self.assertAlmostEqual(0.0, calc_bspline(3, knots, 2, x))
        self.assertAlmostEqual(0.0, calc_deriv_bspline(3, knots, 2, x))

    def test_knots_index_list(self):
        knots = make_knots((0.0, 3), 1.0, (2.5, 2))
        index_list = knots_index_list(knots)
        self.assertEqual(6, len(index_list))
        self.assertEqual(0, index_list[0])
        self.assertEqual(2, index_list[-1])

    def test_non0_q_list(self):
        non0 = non0_q_list(3, self.knots, 0)
        self.assertEqual(4, len(non0))
        self.assertTrue(non0[0])
        for x in non0[1:]:
            self.assertFalse(x)
        non0 = non0_q_list(3, self.knots, 2)
        self.assertEqual([True, True, True, False], non0)

    def test_has_non0(self):
        order = 3
        knots = lin_knots(0.0, 5.0, 6)(order)
        x = 0.5

        def f(i, x):
            return calc_bspline(order, knots, i, x)

        def check(i, j, x):
            self.assertEqual(f(i, x) * f(j, x) != 0.0, has_non0(order, i, j))

        check(0, 1, 0.5)
        check(0, 2, 0.5)
        check(0, 3, 0.5)
        check(2, 5, 1.5)
        check(1, 5, 2.5)

    def test_BSplineSet(self):
        knots = make_knots((0.0, 3), 1.0, 2.0, 3.0, 4.0, (5.0, 3))
        bspline_set = BSplineSet(3, knots)
        self.assertEqual(5, len(bspline_set.basis))

    def test_smat(self):
        bspline_set = BSplineSet(3, lin_knots(0.0, 10.0, 11))
        smat0 = bspline_set.s_mat_old()
        smat1 = bspline_set.s_mat()
        eps = max(flatten((smat0 - smat1).toarray().tolist()))
        self.assertTrue(abs(eps) < 10.0**(-10))

    def test_eri(self):
        compare_flag = True

        bspline_set = BSplineSet(5, lin_knots(0.0, 10.0, 24))
        L = 2
        t0 = time.clock()
        eri_cpp = bspline_set.eri_mat_cpp(L)
        t_cpp = time.clock()-t0

        if compare_flag:
            t0 = time.clock()
            eri_slice = bspline_set.eri_mat(2)
            t_slice = time.clock()-t0

            t0 = time.clock()
            eri_py = bspline_set.eri_mat_old(2)
            t_py = time.clock()-t0

            self.assertEqual(len(eri_py.nonzero()[0]),
                             len(eri_cpp.nonzero()[0]))
            eps = np.max(np.abs((eri_cpp - eri_py).toarray()))
            self.assertTrue(abs(eps) < 10.0**(-10))

            print "test_eri: c++, py, slice"
            print (t_cpp, t_py, t_slice)
        else:
            print t_cpp

    def test_hydrogen_atom(self):
        basis_set = BSplineSet(9, lin_knots(0.0, 100.0, 101))
        d2_mat = basis_set.d2_mat()
        v_mat = basis_set.v_mat(lambda r: -1.0/r)
        s_mat = basis_set.s_mat()
        h_mat = -0.5*d2_mat + v_mat
        (eigs, vecs) = la.eigsh(h_mat, 6, s_mat, sigma=-0.5)
        for (e, n) in zip(eigs, range(1, 6)):
            self.assertAlmostEqual(e, -0.5/(n*n))

if __name__ == '__main__':
    unittest.main()
