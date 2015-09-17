import unittest
from utils import *


class TestUtils(unittest.TestCase):
    def setUp(self):
        pass

    def test_flatten(self):
        self.assertEqual([1, 2, 3, 4], flatten([[1], [2, 3], [4]]))

    def test_uniq(self):
        self.assertEqual([1, 2, 3], uniq([1, 2, 2, 1, 3, 3]))

if __name__ == '__main__':
    unittest.main()
