import numpy as np
from numpy import pi, sqrt
from sympy.physics.wigner import gaunt
from scipy.sparse import bmat
from bspline_bind import ra_inv
from utils import flatten, uniq


def non0_ls(l1, l2):
    """ gives list of qunatum number L which give non 0 <l1,0|l,0|l2,0>
    <l1,m1|l,m|l2,m2> have non 0 value if
    1. l1, l, l2 form triangle
    2. m1 = m + m2
    3. l1 + l + l2 is even integer

    notice that <l1,m1|l,m|l2,m2> is proportional to two 3-j symbol.
    <l1,m1|l,m|l2,m2> prop (L1,L,L2,0,0,0)(L1,L,L2,-M1,M,M2)
    feature 3 is from first 3-j symbol

    Parameters
    ----------
    l1, l2 : non negative Int

    Returns
    -------
    : list of non negative Int

    non0_ls(0, 3)
    >>> [3]
    non0_ls(2, 3)
    >>> [1, 3, 5]
    non0_ls(3, 3)
    >>> [0, 2, 4, 6]
    """
    return [l for l in range(abs(l1-l2), abs(l1+l2)+1) if (l+l1+l2) % 2 == 0]


def mat_h2_plus(bond_length, bspline_set, l_list):
    """ gives hamiltonian matrix and overlap matrix of hydrogem molecule ion

    Parameters
    ----------
    bond_length : Double
         bond length of hydrogen molecular ion
    bspline_set : BSplineSet
    l_list: list of non negative integer
         list of angular quauntum number to use

    Returns
    -------
    h_mat : numpy.ndarray
        hamiltonian matrix
    s_mat : numpy.ndarray
        overlap pmatrix
    """

    # compute matrix in 1D space (B_i|O|B_j)
    # (O=-1/2 d^2/dr^2, 1, 1/r^2, {s^L/g^{L+1} | L<-l_list})
    rs = bspline_set.xs
    d2_mat = bspline_set.d2_mat()
    r2_mat = bspline_set.v_mat(1.0/(rs*rs))
    s_mat = bspline_set.s_mat()
    tmp_L_list = uniq(flatten([non0_ls(L1, L2)
                               for L1 in l_list for L2 in l_list]))
    ra_s_L = {L: bspline_set.v_mat(ra_inv(rs, L, bond_length/2.0))
              for L in tmp_L_list}

    # compute matrix in 3D space (B_iY_L1M1|O|B_jY_L2M2)
    def one_block(L1, L2):
        v = -2.0*sum([sqrt(4.0*pi/(2*L+1)) *
                      float(gaunt(L1, L, L2, 0, 0, 0, 12)) *
                      ra_s_L[L]
                      for L in non0_ls(L1, L2)])
        if L1 == L2:
            L = L1
            t = -0.5 * d2_mat + L*(L+1)*0.5*r2_mat
            return t+v
        else:
            return v

    H_mat = bmat([[one_block(L1, L2)
                   for L1 in l_list]
                  for L2 in l_list])
    S_mat = bmat([[s_mat if L1 == L2 else None
                   for L1 in l_list]
                  for L2 in l_list])
    return (H_mat, S_mat)
