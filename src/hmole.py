import numpy as np
from numpy import pi, sqrt
from sympy.physics.wigner import gaunt
from scipy.sparse import bmat
from bspline_bind import ra_inv
from utils import flatten, uniq
from angmoment import ls_non_zero_YYY


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

    # compute r1 matrix  (B_i|O|B_j)
    # (-1/2 d^2/dr^2, 1, 1/r^2, {s^L/g^{L+1} | L<-l_list})
    rs = bspline_set.xs
    d2_mat = bspline_set.d2_mat()
    r2_mat = bspline_set.v_mat(1.0/(rs*rs))
    s_mat = bspline_set.s_mat()
    tmp_L_list = uniq(flatten([ls_non_zero_YYY(L1, L2)
                               for L1 in l_list for L2 in l_list]))
    #    ra_s_L = {L: bspline_set.v_mat(ra_inv(rs, L, bond_length/2.0))
    #              for L in tmp_L_list}
    ra_s_L = {L: bspline_set.en_mat(L, bond_length/2.0) for L in tmp_L_list}

    # compute r1Y matrix (B_iY_L1M1|O|B_jY_L2M2)
    def one_block(L1, L2):
        v = -2.0*sum([sqrt(4.0*pi/(2*L+1)) *
                      float(gaunt(L1, L, L2, 0, 0, 0, 12)) *
                      ra_s_L[L]
                      for L in ls_non_zero_YYY(L1, L2)])
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
