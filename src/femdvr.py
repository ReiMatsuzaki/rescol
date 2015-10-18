from sympy import Symbol, diff, legendre
from sympy.mpmath import taylor, polyroots
import scipy.sparse as ma
import scipy.sparse.linalg as la
import numpy as np
from utils import with_index, diag_bmat, flatten, repeat, replace_at


"""
fmat: matrix represented by Lobatto Shape basis
chimat: matrix represented by FEM-DVR basis which is continued function
        build from Lobatto shape basis functions
"""


def lin_knots(a, b, num):
    d = (b-a)/(num-1)
    return [a+i*d for i in range(num)]


def exp_knots(a, b, n, gamma):
    g = float(gamma)
    return [a+(b-a)*(np.exp((g*i)/(n-1))-1)/(np.exp(g)-1) for i in range(n)]


def lobgauss(n, a=None, b=None):
    """gives Gauss-Lobatto quadrature points and weights

    Parameters
    ----------
    n : positive integer
      number of quadrature points and weights
    a, b : real value
      left and right edge of region. None means a=-1 and b=+1

    Returns
    -------
    (xs, ws) : ([Double], [Double])
      quadrature points and weights

    Reference
    ---------
    calculation of zero points for polynomial
    http://docs.sympy.org/latest/modules/mpmath/functions/orthogonal.html
    """
    x = Symbol("x")
    lp = diff(legendre(n-1, x), x)
    xs_in = ([-1.0] +
             polyroots(taylor(lambda x0: lp.subs(x, x0), 0, n-2)[::-1]) +
             [+1.0])
    ws_in = [2.0/(n * (n-1) * legendre(n-1, x0) ** 2) for x0 in xs_in]

    if a is None:
        a = -1.0
    if b is None:
        b = 1.0

    xs = [(b-a)/2.0*x+(b+a)/2.0 for x in xs_in]
    ws = [(b-a)/2.0*w for w in ws_in]

    def to_np(ys):
        return np.array(map(float, ys))

    return (to_np(xs), to_np(ws))


def deriv_f(xs, ws, i, m, mp):
    """gives derivative of Lobatto shape function f_{i,m}
    at quadrature point x_{i,mp}"""
    (ne, n) = xs.shape
    if m == mp:
        if m == n-1:
            return +1/(2.0*ws[i, m])
        elif m == 0:
            return -1/(2.0*ws[i, m])
        else:
            return 0
    else:
        ks = range(n)
        ks.remove(m)
        ks.remove(mp)
        return 1.0 / (xs[i, m]-xs[i, mp]) * \
            reduce(lambda a, b: a*b,
                   [(xs[i, mp]-xs[i, k])/(xs[i, m]-xs[i, k]) for k in ks])


def trans_mat_f_to_chi(n, ne):
    zero = ma.coo_matrix([repeat(0, n-2)])
    middle = ma.identity(n-2)
    comb = ma.coo_matrix([[1], [1]])
    ds = flatten(repeat([middle, comb], ne-1)) + [middle]
    nones = repeat(None, 2*(ne-1))
    bmat_args = ([[zero] + nones] +
                 diag_off_none(ds) +
                 [nones + [zero]])
    return ma.bmat(bmat_args)


def diag_off_none(xs):
    n = len(xs)
    return [replace_at(repeat(None, n), i, x)
            for (x, i) in with_index(xs)]


def f_to_chi(m, n, ne):
    t = trans_mat_f_to_chi(n, ne)
    return t.T.dot(m.dot(t))


def lobatto_shape(xs, m, x):
    """gives value of m th LobattoShape function for quadrature points xs"""
    if x < xs[0]:
        return 0.0
    if xs[-1] < x:
        return 0.0
    xs_nonm = [xi for (xi, mi) in with_index(xs) if mi != m]
    return reduce(lambda a, b: a*b, [(x-xi)/(xs[m]-xi) for xi in xs_nonm])


class FemDvrSet:
    """ FEM-DVR set"""
    def __init__(self, n, ts):
        self.n = n
        self.ts = ts
        xws_list = np.array([lobgauss(n, a, b)
                             for (a, b)
                             in zip(ts[:-1], ts[1:])])

        # by function 'rollaxis', the shape is changed to
        # xws_list.shape = (len(ts)-1, 2, n) -> (2, len(ts)-1, n)
        [self.xs, self.ws] = np.rollaxis(xws_list, 1)

    def prepare_deriv(self):
        # compute derivative of Lobatto-Shapiro functions
        # fp_list.shape = (# of element, # of basis in ele, # of quad point)
        ne = len(self.ts)-1
        self.fp_list = np.array(
            [[[deriv_f(self.xs, self.ws, i, m, mp)
               for mp in range(self.n)]
              for m in range(self.n)]
             for i in range(ne)]
        )

    def s_fmat(self):
        # diags : gives diagonal matrix
        # ravel : flatten numpy array
        return ma.diags([np.ravel(self.ws)], [0])

    def v_fmat(self, v):
        vs = np.r_[np.array([0.0]), v(np.ravel(self.xs)[1:])]
        return ma.diags([np.ravel(self.ws) * vs], [0])

    def d2_fmat(self):
        self.prepare_deriv()
        bmat_list = [ma.coo_matrix([[
            -np.sum(fpi*fpj*ws)
            for fpi in fps] for fpj in fps])
                     for (ws, fps) in zip(self.ws, self.fp_list)]
        return diag_bmat(bmat_list)

    def f_to_chi(self, m):
        return f_to_chi(m, self.n, len(self.ts)-1)

    def basis_psi(self, xs):
        return [[lobatto_shape(quad_xs, m, x)
                 for x in xs]
                for quad_xs in self.xs
                for m in range(len(quad_xs))]
