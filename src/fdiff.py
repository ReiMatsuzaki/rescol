import numpy as np
from scipy.sparse import coo_matrix
from utils import flatten
from bspline_bind import ra_inv


class FDiffSet:
    def __init__(self, num, x_max):
        self.h = x_max/num
        self.xs = np.array([self.h*(i+1) for i in range(num)])

    def v_mat(self, v_x):
        v = v_x(self.xs)
        col = np.array(range(len(v)))
        row = np.array(range(len(v)))
        return coo_matrix((v, (row, col)))

    def d2_mat(self):
        num = len(self.xs)
        h2 = self.h**2
        [row, col, dat] = np.array(filter(
            lambda x: 0 <= x[0] and x[0] < num,
            flatten([[(n-1, n, +1/h2),
                      (n, n, -2/h2),
                      (n+1, n, +1/h2)] for n in range(num)]))).T
        return coo_matrix((dat, (row, col)))

    def en_mat(self, L, a):
        return self.v_mat(lambda rs: ra_inv(rs, L, a))

    def eri_mat(self, L):
        print "Not implemented yet."
        raise
