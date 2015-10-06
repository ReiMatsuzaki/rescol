import unittest
import os
import scipy.linalg as la
import numpy as np

class TestDiag(unittest.TestCase):
    def setUp(self):
        pass

    def test_diag(self):
        # compute reference values
        H = [[1.1, 1.2], [1.2, 1.3]]
        S = [[1.0, 0.0], [0.0, 2.0]]
        (vals, vecs) = la.eigh(H, S)

        # computation and get results
        os.system('./diag.out')
        with open('diag_eigvals.dat') as f:
            f.readline() # skip header
            line_list = f.readlines()
        vals_list = [map(float, line.split(" "))[0]
                    for line in line_list]

        # sort
        vals_ref = sorted(vals)
        vals_calc= sorted(vals_list)

        # test
        for (ref, calc) in zip(vals_ref, vals_calc):
            self.assertTrue(
                abs(ref-calc)<0.00001,
            "{0}, {1}".format(ref, calc))

if __name__ == '__main__':
    unittest.main()

