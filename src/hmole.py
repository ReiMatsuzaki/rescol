import numpy as np
from numpy import pi, sqrt
from scipy.sparse import bmat, coo_matrix
from utils import flatten, uniq, synthesis_mat
from angmoment import ls_non_zero_YYY, y1mat_Yqk


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
    d2_rmat = bspline_set.d2_mat()
    r2_rmat = bspline_set.v_mat(1.0/(rs*rs))
    s_rmat = bspline_set.s_mat()
    tmp_L_list = uniq(flatten([ls_non_zero_YYY(L1, L2)
                               for L1 in l_list for L2 in l_list]))
    en_r1mat_L = {L: bspline_set.en_mat(L, bond_length/2.0)
                  for L in tmp_L_list}

    # compute y1 matrix (Y_L1|P_L(w_A)|Y_L2)
    en_y1mat_L = {L: coo_matrix([[np.sqrt(4.0*np.pi/(2*L+1)) *
                                  y1mat_Yqk((L1, 0), (L, 0), (L2, 0))
                                  for L1 in l_list]
                                 for L2 in l_list])
                  for L in tmp_L_list}
    LL_y1mat = coo_matrix(np.diag([1.0*L*(L+1) for L in l_list]))
    diag_y1mat = coo_matrix(np.diag([1 for L in l_list]))

    # compute r1y1 matrix
    h_r1y1mat = (synthesis_mat(-0.5*d2_rmat, diag_y1mat) +
                 synthesis_mat(+0.5*r2_rmat, LL_y1mat) - 2.0 *
                 sum([synthesis_mat(en_r1mat_L[L], en_y1mat_L[L])
                      for L in tmp_L_list]))
    s_r1y1mat = synthesis_mat(s_rmat, diag_y1mat)

    return (h_r1y1mat, s_r1y1mat)


def mat_h2_plus_old(bond_length, bspline_set, l_list):
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
    d2_rmat = bspline_set.d2_mat()
    r2_rmat = bspline_set.v_mat(1.0/(rs*rs))
    s_rmat = bspline_set.s_mat()
    tmp_L_list = uniq(flatten([ls_non_zero_YYY(L1, L2)
                               for L1 in l_list for L2 in l_list]))
    en_r1mat_L = {L: bspline_set.en_mat(L, bond_length/2.0)
                  for L in tmp_L_list}

    # compute r1Y matrix (B_iY_L1M1|O|B_jY_L2M2)
    def one_block(L1, L2):
        v = -2.0*sum([sqrt(4.0*pi/(2*L+1)) *
                      y1mat_Yqk((L1, 0), (L, 0), (L2, 0)) *
                      en_r1mat_L[L]
                      for L in ls_non_zero_YYY(L1, L2)])
        if L1 == L2:
            L = L1
            t = -0.5 * d2_rmat + L*(L+1)*0.5*r2_rmat
            return t+v
        else:
            return v

    H_mat = bmat([[one_block(L1, L2)
                   for L1 in l_list]
                  for L2 in l_list])
    S_mat = bmat([[s_rmat if L1 == L2 else None
                   for L1 in l_list]
                  for L2 in l_list])
    return (H_mat, S_mat)
