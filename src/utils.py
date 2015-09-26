from scipy.sparse import coo_matrix


def flatten(xs):
    return reduce(lambda a, b: a+b, xs)


def uniq(xs):
    return list(set(xs))


def with_index(xs):
    return zip(xs, range(len(xs)))


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

    col = [i+na_col*j for j in b.col for i in a.col]
    row = [k+na_row*l for l in b.row for k in a.row]
    data = [ik*jl for jl in b.data for ik in a.data]
    return coo_matrix((data, (row, col)))
