import unittest
from SignAlgebra.AMinus import *


class TestAMinus(unittest.TestCase):

    def test_mul(self):
        am = AMinus([-1, 1, 1, -1])
        r = am.polyring
        U1 = am.polyring['U1']
        U2 = am.polyring['U2']

        # First multiplication of single elements
        # double crosses black, should be empty product
        elt1 = am.generator({4: 0, 3: 1})
        elt2 = am.generator({0: 4, 1: 3})
        assert elt1 * elt2 == am.zero()

        # double crosses two rightward strands, should be nontrivial product
        elt3 = am.generator({4: 0, 3: 1})
        elt4 = am.generator({0: 0, 1: 3})
        assert elt3 * elt4 == (U1 * U2) * am.generator({4: 0, 3: 3})

        # just straightforward multiplication
        elt5 = am.generator({4: 1, 3: 2})
        elt6 = am.generator({2: 1, 1: 0})
        assert elt5 * elt6 == am.generator({4: 0, 3: 1})

    def test_diff(self):
        am = AMinus([-1, 1, 1, -1])
        r = am.polyring
        U1 = am.polyring['U1']
        U2 = am.polyring['U2']

        elt = am.generator({0: 4, 1: 2, 3: 1})
        assert elt.diff() == U1 * am.generator({0: 4, 1: 1, 3: 2}) + am.generator({0: 2, 1: 4, 3: 1})

    def test_gradings(self):
        pass  # TODO
