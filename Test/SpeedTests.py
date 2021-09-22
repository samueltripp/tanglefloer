from Modules.CTMinus import *
from Modules.Module import *
from Modules.ETangleStrands import *
import unittest
import timeit

from Tangles import TangleRenderer


class TestSpeed(unittest.TestCase):
    # - 130/10
    def test_da_speed(self):
        cup = ETangle(ETangle.Type.OVER, (-1, 1, 1), 1)
        for _ in range(10):
            da = type_da(cup)
        self.assertTrue(True)


if __name__ == '__main__':
    unittest.main()
