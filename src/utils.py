from scipy.sparse import coo_matrix
import numpy as np
import bspline_bind


def flatten(xs):
    return reduce(lambda a, b: a+b, xs)


def uniq(xs):
    return list(set(xs))


def with_index(xs):
    return zip(xs, range(len(xs)))


def outer_sum(a, b):
    """ gives summation version of outer product

    {a_i+b_j}_{i, j}

    Parameters
    ----------
    a, b : 1D numpay.array

    Results
    -------
    out : 1D numarray.array
    """
    return bspline_bind.outer_sum(a, b)


def synthesis_mat_old(a, b):
    """ gives conbination of two matrix a and b

    Parameters
    ----------
    a, b : scipy.sparse.occ_matrix
       a and b must be square matrix and

    Returns
    -------
    out : scipy.sparse.occ_matrix
    """
    (na_row, na_col) = a.shape
    (nb_row, nb_col) = b.shape

    n_row = na_row * nb_row
    n_col = na_col * nb_col

    col = [i+na_col*j for j in b.col for i in a.col]
    row = [k+na_row*l for l in b.row for k in a.row]
    data = [ik*jl for jl in b.data for ik in a.data]
    return coo_matrix((data, (row, col)))


def synthesis_mat(a, b):
    """ gives conbination of two matrix a and b

    Parameters
    ----------
    a, b : scipy.sparse.occ_matrix
       a and b must be square matrix and

    Returns
    -------
    out : scipy.sparse.occ_matrix
    """
    (na_row, na_col) = a.shape
    (nb_row, nb_col) = b.shape

    n_row = na_row * nb_row
    n_col = na_col * nb_col

    num = len(a.data)*len(b.data)

    ncol = na_col*b.col
    col = bspline_bind.outer_sum(ncol, a.col)
    nrow = na_row*b.row
    row = outer_sum(nrow, a.row)

    data = np.reshape(np.outer(b.data, a.data), (num))
    return coo_matrix((data, (row, col)), shape=(n_row, n_col))


def write_coo_mat(mat, file_name):
    mat = coo_matrix(mat)
    row = mat.row
    col = mat.col
    dat = mat.data
    with open(file_name, mode='w') as f:
        for (r, c, d) in zip(row, col, dat):
            f.write("{0} {1} {2}\n".format(r, c, d))


def read_coo_mat(file_name):
    [row, col, dat] = np.loadtxt(file_name,
                                 delimiter=" ",
                                 usecols=(0, 1, 2)).T
    return coo_matrix((dat, (row, col)))
