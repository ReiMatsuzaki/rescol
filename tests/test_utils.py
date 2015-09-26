import unittest
from utils import *
import numpy as np
from scipy.sparse import coo_matrix


class TestUtils(unittest.TestCase):
    def setUp(self):
        pass

    def test_flatten(self):
        self.assertEqual([1, 2, 3, 4],
                         flatten([[1], [2, 3], [4]]))

    def test_uniq(self):
        self.assertEqual([1, 2, 3],
                         uniq([1, 2, 2, 1, 3, 3]))

    def test_synthesis_mat(self):
        """
        m1 =
        [1, 0, 2, -1]
        [0, 0, 3,  0]
        [4, 5, 6, -2]
        m2 =
        [1,-2, -5]
        [0, 3, -6]
        """
        row = np.array([0, 0, +0, 1, 2, 2, 2, 2])
        col = np.array([0, 2, +3, 2, 0, 1, 2, 3])
        dat = np.array([1, 2, -1, 3, 4, 5, 6, -2])
        m1 = coo_matrix((dat, (row, col)))
        row = np.array([+0, +0, +0, +1, +1])
        col = np.array([+0, +1, +2, +1, +2])
        dat = np.array([+1, -2, -5, +3, -6])
        m2 = coo_matrix((dat, (row, col)))

        self.assertEqual((3, 4), m1.shape)
        self.assertEqual((2, 3), m2.shape)

        m3 = synthesis_mat(m1, m2).toarray()
        self.assertEqual((6, 12), m3.shape)
        self.assertAlmostEqual(2.0, m3[0, 2])
        self.assertAlmostEqual(3.0, m3[1, 2])
        self.assertAlmostEqual(-6.0, m3[1, 4+2])
        self.assertAlmostEqual(+9.0, m3[3+1, 4+2])


if __name__ == '__main__':
    unittest.main()
