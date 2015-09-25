import numpy as np
import sympy.physics.wigner as wigner

"""
Yqk,Y_LM : Spherical Harmonics
y1mat    : matrix or matrix element of Y_LM
y1redmat : reduced matrix element or its element of Y_LM
coupledY : combination of two Y_LM
y2mat    : matrix or its element of coupledY
y2redmat : reduced matrix or its element of coupledY
Pq_r12   : P_q(cos r_12) = 4pi/(2q+1) (-)^k Y_q-kY_qk
"""


def prod_LLL(*args):
    """ gives product of (2 L_i+1).

    Parameters
    ----------
    *args: array_like of number

    Returns
    -------
    out : number
    """
    return reduce(lambda a, b: a*b, [2*L+1 for L in args])


def N3j(*args):
    """ gives numerical value of 3-j symbol
    """
    return wigner.wigner_3j(*args)


def N6j(*args):
    """ gives numerical value of 6-j symbol
    """
    return wigner.wigner_6j(*args)


def triangle_q(L1, L2, L3):
    """ gives whether given three numbers form triangle or not.
    """
    return (
        (abs(L1-L2) <= L3 and L3 <= (L1+L2)) or
        (abs(L2-L3) <= L1 and L1 <= (L2+L3)) or
        (abs(L3-L1) <= L2 and L2 <= (L3+L1)))


def zero_YYY_q(L1, L2, L3):
    """ gives <Y(L1)||Y(L2)||Y(L3)> == 0
    """
    return (not triangle_q(L1, L2, L3) or
            (L1+L2+L3) % 2 != 0)


def _y1redmat_Yq(L1, q, L2):
    return (np.sqrt(prod_LLL(L1, q, L2)/(4*np.pi)) *
            N3j(L1, q, L2, 0, 0, 0))


def y1redmat_Yq(L1, q, L2):
    """ gives <Y(L1)||Y(q)||Y(L2)>
    """
    if zero_YYY_q(L1, q, L2):
        return 0
    else:
        return _y1redmat_Yq(L1, q, L2)


def _y1mat_Yqk(L1M1, qk, L2M2):
    (L1, M1) = L1M1
    (q, k) = qk
    (L2, M2) = L2M2
    return ((-1)**(L1-M1) *
            y1redmat_Yq(L1, q, L2) *
            N3j(L1, q, L2,
                -M1, k, M2))


def y1mat_Yqk(L1M1, qk, L2M2):
    """ gives <Y_L1M1|Y_qk|YL2M2>
    """
    (L1, M1) = L1M1
    (q, k) = qk
    (L2, M2) = L2M2
    if zero_YYY_q(L1, q, L2):
        return 0
    if M1 != k + M2:
        return 0
    return _y1mat_Yqk(L1M1, qk, L2M2)


class CoupledY:
    def __init__(self, L1L2, L, M):
        (L1, L2) = L1L2
        if(not triangle_q(L1, L2, L)):
            msg = """
            (L1, L2, L) must form triangle relations."
            (L1, L2, L) = ({0}, {1}, {2})
            """.format(L1, L2, L)
            raise Exception(msg)
        if(abs(M) > L):
            msg = """
            |M|>L detected. \n(L,M)=({0},{1})"
            """.format(L, M)
            raise Exception(msg)
        self.L1 = L1
        self.L2 = L2
        self.L = L
        self.M = M
        print self.L

    def __str__(self):
        return """
        Y[({0},{1}) {2},{3}]
        """.format(self.L1, self.L2, self.L, self.M)

    def __repr__(self):
        return self.__str__()


def _y2mat_YYq(yp, q, y):
    L = y.L
    M = y.M
    t1 = (-1)**(yp.L2+L+y.L1)
    t2 = y1redmat_Yq(yp.L2, q, y.L2)
    t3 = y1redmat_Yq(yp.L1, q, y.L1)
    t4 = N6j(yp.L1, yp.L2, L,
             y.L2,  y.L1,  q)
    return t1*t2*t3*t4


def y2mat_YYq(yp, q, y):
    """ gives matrix element of scalar tensor YY
    YY = sum_k (-)^qY_q-kY_k

    Parameters
    ----------
    yp, y : CoupledY
    q     : non negative integer

    Results
    -------
    out : real number
    """
    if ((yp.L != y.L or yp.M != y.M or
         zero_YYY_q(yp.L1, q, y.L1) or
         zero_YYY_q(yp.L2, q, y.L2))):
        return 0
    return _y2mat_YYq(yp, q, y)


def y2mat_Pq_r12(yp, q, y):
    """ gives matrix element of P_q(cos theta_12)

    P_q: Legendre function of order q
    theta_12 : angle between two electron

    P_q = 4pi/(2q+1) YY
    """
    return 4.0*np.pi/(2*q+1)*y2mat_YYq(yp, q, y)


def _y2redmat_Y1q(yp, q, y):
    if yp.L2 != y.L2:
        return 0
    t1 = (-1)**(y.L+yp.L1+y.L2+q)
    t2 = np.sqrt(prod_LLL(y.L, yp.L))
    t3 = y1redmat_YL(yp.L1, q, y.L1)
    t4 = N6j(yp.L1, yp.L, y.L2,
             y.L,   y.L1, q)
    return t1*t2*t3*t4


def y2redmat_Y1q(yp, q, y):
    """ gives reduced matrix element of Y_q(1)

    Parameters
    ----------
    yp, y : CoupledY
    q : non negative integer

    Results
    -------
    out : real number
    """
    if ((yp.L2 != y.L2 or
         zero_YYY_q(yp.L1, q, y.L1))):
        return 0
    return _y2redmat_Y1q(yp, q, y)


def _y2mat_Y1qk(yp, q, k, y):
    t = (-1)**(yp.L-yp.M)
    rm = y2redmat_Y1q(yp, q, y)
    tj = N3j(+yp.L, +q, +y.L,
             -yp.M, +k, +y.M)
    return t * rm * tj


def y2mat_Y1qk(yp, q, k, y):
    """ gives matrix elemtns of Y_qk(1)
    """
    if ((y.L2 != yp.L2 or
         zero_YYY_q(yp.L1, q, y.L1) or
         not triangle_q(yp.L, q, y.L) or
         yp.M != k + y.M)):
        return 0
    t = (-1)**(yp.L-yp.M)
    rm = y2redmat_Y1q(yp, q, y)
    tj = N3j(+yp.L, +q, +y.L,
             -yp.M, +k, +y.M)
    return t * rm * tj


def y2mat_Pq_r1A(yp, q, y):
    """ gives matrix elemtn of P_q(cos theta_rA)

    A is on positive z-axis
    """
    return np.sqrt(4.0*np.pi/(2*q+1)) * y2mat_Y1qk(yp, q, 0, y)
