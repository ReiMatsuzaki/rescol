from sympy import Symbol, diff, legendre
from sympy.mpmath import taylor, polyroots
import scipy.sparse as ma
import scipy.sparse.linalg as la
import numpy as np
from utils import with_index


"""
fmat: matrix represented by Lobatto Shape basis
chimat: matrix represented by FEM-DVR basis which is continued function
        build from Lobatto shape basis functions
"""


def lin_knots(a, b, num):
    d = (b-a)/(num-1)
    return [a+i*d for i in range(num)]


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


def f_to_chi(m, n, ne):
    pass


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
        # fp_list.shape = ()
        ne = len(ts)-1
        self.fp_list = np.array(
            [[[deriv_f(self.xs, self.ws, i, m, mp)
               for mp in range(n)]
              for m in range(n)]
             for i in range(ne)])

    def s_fmat(self):
        # diags : gives diagonal matrix
        # ravel : flatten numpy array
        return ma.diags([np.ravel(self.ws)], [0])

    def v_fmat(self, v):
        vs = v(np.ravel(self.xs))
        return ma.diags([np.ravel(self.ws) * vs], [0])

    def d2_fmat(self, v):
        self.prepare_deriv()

        bmat_list = [ma.coo_matrix([[
            -np.sum(fpi*fpj*self.ws)
            for fpi in fps] for fpj in fps])
                for fps in self.fp_list]
        return diag_bmat(bmat_list)
