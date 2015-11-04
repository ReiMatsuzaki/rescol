import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import coo_matrix
from bspline_bind import calc_bspline_xs, calc_deriv_bspline_xs
from bspline_bind import ra_inv, dot_abwAcdw, eri_mat
from utils import *


def make_knots(*args):
    """ from arguments, make knots data
    knots data is
    [(KnotLocation, Number)]:: [(Double, Int)]
    args are (Doube, Int) or Double
    make_knots(2.0, 3.0)
    >>> [(2.0,1), (3.0,1)]

    make_knots((2.0,3), 1.5)
    >>> [(2.0,3), (1.5,1)]
    """
    return [x if type(x) == tuple else (x, 1) for x in args]


def lin_knots(x0, x1, num):

    def func(order):
        xs = [x0 + (x1-x0)/(num-1)*i for i in range(num)]
        xn_list = [(xs[0], order)] + xs[1:-1] + [(xs[-1], order)]
        return make_knots(*xn_list)
    return func


def exp_knots(x0, x1, num, gamma):
    a = float(x0)
    b = float(x1)
    g = float(gamma)
    n = num
    def func(order):
        xs = [a+(b-a)*(np.exp((g*i)/(n-1))-1)/(np.exp(g)-1) for i in range(n)]
        xn_list = [(xs[0], order)] + xs[1:-1] +[(xs[-1], order)]
        return make_knots(*xn_list)
    return func

def make_intervals(knots):
    """ from knots, make intervals of knots.

    intervals = [(Double,Double)]

    make_intervals([(1.1,3), (1.2,2), (3.1,1)])
    >>> [ (1.1, 1.2), (1.2,3.1) ]
    """
    xs = [x for (x, n) in knots]
    return [(a, b) for (a, b) in zip(xs[0:-1], xs[1:])]


def make_ts(knots):
    """ gives {t_i}_i from knots
    make_ts([(1.1,3), (2.0,1), (3.4,3)])
    >>> array([1.1, 1.1, 1.1, 2.0, 3.4, 3.4, 3.4])
    """
    return np.array(reduce(lambda a, b: a+b, [[x]*n for (x, n) in knots]))


def calc_bspline(order, knots, index, x):
    """ compute value of B-Spline function at x
    order(Int) : B-Spline order (k in paper)
    knots( [(double,Int)]) : knot location and multiplicity
    index(Int) : index of B-Spline (0 origin)
    x(double) :
    """
    ts = make_ts(knots)
    xs = np.array([x])
    ys = calc_bspline_xs(order, ts, index, xs)
    return ys[0]


def calc_deriv_bspline(order, knots, index, x):
    """ compute derivative value of B-Spline function at x
    order(Int) : B-Spline order (k in paper)
    knots( [(double,Int)]) : knot location and multiplicity
    index(Int) : index of B-Spline (0 origin)
    x(double) :
    """
    ts = make_ts(knots)
    xs = np.array([x])
    ys = calc_deriv_bspline_xs(order, ts, index, xs)
    return ys[0]


def knots_index_list(knots):
    """ return index list of knots point
    knots_index_list( [(0.0, 3), (1.0, 1), (2.0,2)] )
    >>> [0,0,0,1,2,2]
    """
    def _knots_index_list(res, index, knots):
        if len(knots) == 0:
            return res
        (x, num) = knots[0]
        res1 = [index for i in range(num)]
        return _knots_index_list(res + res1,
                                 index + 1,
                                 knots[1:])
    return _knots_index_list([], 0, knots)


def non0_q_list(order, knots, index):
    """ function have non 0 value or not list

    knots = make_knots((0.0,3), 1.0, 2.0, 3.0, 4.0)
    non0_q_list(3, knots, 0)
    >>> [true, false, false, false]
    non0_q_list(3, knots, 1)
    >>> [true, true, false, false]
    non0_q_list(3, knots, 3)
    >>> [false, false, true, true]
    """
    index_list = knots_index_list(knots)
    i0 = index_list[index]
    i1 = index_list[index+order]
    return [False]*i0 + [True]*(i1-i0) + [False]*(index_list[-1]-i1)


def has_non0(order, i, j):
    """ multiplicative of two B-Spline functuon B_i(x)B_j(x) have
    non 0 value or not.
    """
    return abs(i-j) < order


def num_bspline(order, knots):
    """ return number of correct B-Spline functions
    """
    return len(knots)-1+(order-1)


class BSpline:
    """ B-Spline basis function
    index: Int    index of this B-Spline function
    non0: [Bool]  true=>have non 0 value in that interval
    val: [Double] list of function value with weight factor
    deriv: [Double] list of derivative of function value with weight factor
    """
    def __init__(self, index, non0_q, val, deriv):
        self.index = index
        self.non0_q = non0_q
        self.val = val
        self.deriv = deriv


class BSplineSet:
    """ B-Spline basis set for compute matrix element,
    vector element, x-representation
    """
    def __init__(self, order, knots):
        """ Constructor for B-Spline set.
        order: Int
        knots: (Int)->[(Double,Int)] or [(Double,Int)]
        """
        self.knots = knots if type(knots) == list else knots(order)
        self.order = order
        intervals = make_intervals(self.knots)
        (xs_quad, ws_quad) = leggauss(order)
        xs = np.array(
            [(b-a)*x/2.0+(b+a)/2.0 for (a, b) in intervals for x in xs_quad])
        ws = np.array(
            [w*(b-a)/2.0 for (a, b) in intervals for w in ws_quad])
        ts = make_ts(self.knots)

        def one_basis(i):
            non0 = non0_q_list(order, self.knots, i)
            val = calc_bspline_xs(order, ts, i, xs)
            deriv = calc_deriv_bspline_xs(order, ts, i, xs)
            return BSpline(i, non0, val, deriv)

        self.xs = xs
        self.ws = ws
        self.ts = ts
        self.basis = [one_basis(i) for i
                      in range(1, num_bspline(self.order, self.knots)-1)]

    def calc_coo_matrix(self, mat_ele):
        """ TO BE REMOVED"""
        return self.calc_coo_matrix_old(mat_ele)

    def calc_mat(self, ab_to_ys):

        def ele(a, b):
            (i0, i1) = self.non0_index(a, b)
            ys_list = ab_to_ys(a, b)
            return np.dot(self.ws[i0:i1],
                          reduce(lambda t, u: t*u,
                                 [ys[i0:i1] for ys in ys_list]))

        (row, col, data) = np.array(
            [[i, j, ele(a, b)]
             for (a, i) in with_index(self.basis)
             for (b, j) in with_index(self.basis)
             if has_non0(self.order, a.index, b.index)]).T
        n = len(self.basis)
        return coo_matrix((data, (row, col)), shape=(n, n))

    def calc_coo_matrix_old(self, mat_ele):
        """TOO BE REMOVED
        Gives COO matrix using B-Spline function
        Parameters
        ----------
        mat_ele : (BSpline,BSpline) -> [Double]
                function which compute discrete representation
        Outpus
        -------
        CSR_matrix
        """
        n = len(self.basis)
        (row, col, data) = np.array(
            [[i, j, np.dot(self.ws, mat_ele(a, b))]
             for (a, i) in zip(self.basis, range(n))
             for (b, j) in zip(self.basis, range(n))
             if has_non0(self.order, a.index, b.index)]).T
        return coo_matrix((data, (row, col)), shape=(n, n))

    def non0_index(self, a, b):
        """give matrix element between BSpline a and b """
        non0_q = [aq and bq for (aq, bq) in zip(a.non0_q, b.non0_q)]
        non0_is = [i for (q, i) in with_index(non0_q) if q]
        i0 = non0_is[0]*self.order
        i1 = (non0_is[-1]+1)*self.order
        return (i0, i1)

    def d2_mat_dense(self):
        return np.array([[-np.sum(self.ws*a.deriv*b.deriv)
                          for a in self.basis] for b in self.basis])

    def d2_mat(self):
        return -self.calc_mat(lambda a, b: [a.deriv, b.deriv])

    def d2_mat_old(self):
        n = len(self.basis)
        (row, col, data) = np.array(
            [[i, j, -np.sum(self.ws*a.deriv*b.deriv)]
             for (a, i) in zip(self.basis, range(n))
             for (b, j) in zip(self.basis, range(n))
             if has_non0(self.order, a.index, b.index)]).T
        return coo_matrix((data, (row, col)), shape=(n, n))

    def s_mat_old(self):
        return self.calc_coo_matrix_old(lambda a, b: a.val*b.val)

    def s_mat(self):
        return self.calc_mat(lambda a, b: [a.val, b.val])

    def v_mat(self, v_x):
        v = v_x if type(v_x) == np.ndarray \
            else np.array([v_x(x) for x in self.xs])
        return self.calc_mat(lambda a, b: [a.val, v, b.val])

    def en_mat(self, L, a):
        return self.v_mat(ra_inv(self.xs, L, a))

    def two_v_mat(self, v):

        def one_ele(a, b, c, d):
            (i0, i1) = self.non0_index(a, c)
            (j0, j1) = self.non0_index(b, d)
            return np.dot(self.ws[i0:i1]*a.val[i0:i1]*c.val[i0:i1],
                          np.dot(v[i0:i1, j0:j1],
                                 self.ws[j0:j1]*b.val[j0:j1]*d.val[j0:j1]))

        usus = [(a, b) for a in self.basis for b in self.basis]
        [row, col, data] = np.array(
            [[I, J, one_ele(a, b, c, d)]
             for ((a, b), I) in with_index(usus)
             for ((c, d), J) in with_index(usus)
             if has_non0(self.order, a.index, c.index) and
             has_non0(self.order, b.index, d.index)]).T
        return coo_matrix((data, (row, col)), shape=(len(usus), len(usus)))

    def two_v_mat_old(self, v):

        def one_ele(a, b, c, d):
            return np.dot(self.ws*a.val*c.val,
                          np.dot(v, self.ws*b.val*d.val))

        usus = [(a, b) for a in self.basis for b in self.basis]
        [row, col, data] = np.array(
            [[I, J, one_ele(a, b, c, d)]
             for ((a, b), I) in with_index(usus)
             for ((c, d), J) in with_index(usus)
             if has_non0(self.order, a.index, c.index) and
             has_non0(self.order, b.index, d.index)]).T
        return coo_matrix((data, (row, col)), shape=(len(usus), len(usus)))

    def eri_mat(self, L):
        return self.eri_mat_cpp(L)

    def eri_mat_cpp(self, L):
        n = len(self.basis)
        bs_vals = np.hstack([u.val for u in self.basis])
        (data, row, col) = eri_mat(bs_vals, self.xs,
                                   self.ws, L, self.order)

        return coo_matrix((data, (row, col)),
                          shape=(n*n, n*n))

    def eri_mat_dense(self, L):

        def eri(x, y, L):
            s = min(x, y)
            g = max(x, y)
            return (s**L)/(g**(L+1))

        eri_ij = np.array([[eri(x, y, L)
                            for x in self.xs]
                           for y in self.xs])

        def one_ele(a, b, c, d):
            vi = self.ws*a.val*c.val
            vj = self.ws*b.val*d.val
            return np.dot(vi, np.dot(eri_ij, vj))
        usus = [(a, b) for a in self.basis for b in self.basis]
        return [[one_ele(a, b, c, d) for (a, b) in usus]
                for (c, d) in usus]

    def m_vec(self, m_x):
        m = np.array([m_x(x) for x in self.xs])
        return np.array([np.sum(self.ws*m*a.val) for a in self.basis])

    def psi(self, cs, xs):
        ys_list = np.array([c*calc_bspline_xs(self.order, self.ts, u.index, xs)
                            for (c, u) in zip(cs, self.basis)])
        return np.sum(ys_list, axis=0)

    def basis_psi(self, xs):
        return np.array([calc_bspline_xs(self.order, self.ts, u.index, xs)
                         for u in self.basis])
