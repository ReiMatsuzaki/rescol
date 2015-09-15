import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.sparse import csr_matrix
from bspline_bind import calc_bspline_xs, calc_deriv_bspline_xs


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
    def __init__(self, index, val, deriv):
        self.index = index
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
            val = calc_bspline_xs(order, ts, i, xs)
            deriv = calc_deriv_bspline_xs(order, ts, i, xs)
            return BSpline(i, val, deriv)

        self.xs = xs
        self.ws = ws
        self.ts = ts
        self.basis = [one_basis(i) for i 
                      in range(1, num_bspline(self.order, self.knots)-1)]

    def calc_csr_matrix(self, mat_ele):
        """Gives CSR matrix using B-Spline function
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
            [[i, j, np.sum(self.ws*mat_ele(a, b))]
             for (a, i) in zip(self.basis, range(n))
             for (b, j) in zip(self.basis, range(n))
             if has_non0(self.order, a.index, b.index)]).T
        return csr_matrix((data, (row, col)), shape=(n, n))

    def d2_mat_dense(self):
        return np.array([[-np.sum(self.ws*a.deriv*b.deriv)
                          for a in self.basis] for b in self.basis])

    def d2_mat(self):
        return self.calc_csr_matrix(lambda a, b: -a.deriv*b.deriv)

    def d2_mat_old(self):
        n = len(self.basis)
        (row, col, data) = np.array(
            [[i, j, -np.sum(self.ws*a.deriv*b.deriv)]
             for (a, i) in zip(self.basis, range(n))
             for (b, j) in zip(self.basis, range(n))
             if has_non0(self.order, a.index, b.index)]).T
        return csr_matrix((data, (row, col)), shape=(n, n))

    def s_mat(self):
        return self.calc_csr_matrix(lambda a, b: a.val*b.val)

    def v_mat(self, v_x):
        if type(v_x) == np.ndarray:
            v = v_x
        else:
            v = np.array([v_x(x) for x in self.xs])
        return self.calc_csr_matrix(lambda a, b: a.val*v*b.val)

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
