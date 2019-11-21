from Tangles.Tangle import *
from Tangles.TangleRenderer import *
import unittest


# some simple examples
cup = ETangle(ETangle.Type.CUP, (1, -1), 1)
over = ETangle(ETangle.Type.OVER, (1, -1), 1)
under = ETangle(ETangle.Type.UNDER, (-1, 1), 1)
cap = ETangle(ETangle.Type.CAP, (1, -1), 1)
unknot = Tangle((cup, cap))
unknot2 = Tangle((cup, over, under, cap))

# the trefoil from the paper
t1 = ETangle(ETangle.Type.CUP, (-1, 1), 1)
t2 = ETangle(ETangle.Type.CUP, (-1, 1, -1, 1), 3)
t3 = ETangle(ETangle.Type.OVER, (-1, 1, -1, 1), 2)
t4 = ETangle(ETangle.Type.UNDER, (-1, -1, 1, 1), 1)
t5 = ETangle(ETangle.Type.OVER, (-1, -1, 1, 1), 2)
t6 = ETangle(ETangle.Type.CAP, (-1, 1, -1, 1), 1)
t7 = ETangle(ETangle.Type.CAP, (-1, 1), 1)
trefoil = Tangle((t1, t2, t3, t4, t5, t6, t7))

# more examples
b1 = ETangle(ETangle.Type.UNDER, (1, 1, 1), 1)
b2 = ETangle(ETangle.Type.UNDER, (1, 1, 1), 2)
twist = b1 + b2
torus_knot = twist + twist + twist + twist


class TestTangles(unittest.TestCase):

    def test_signs(self):
        self.assertEqual(cup.left_signs(), (None,))
        self.assertEqual(cup.middle_signs(), (None, 1, -1))
        self.assertEqual(cup.right_signs(), (None, 1, -1))
        self.assertEqual(cap.left_signs(), (None, 1, -1))
        self.assertEqual(cap.middle_signs(), (None, 1, -1))
        self.assertEqual(cap.right_signs(), (None,))
        self.assertEqual(over.left_signs(), (None, 1, -1))
        self.assertEqual(over.middle_signs(), (None, 1, -1))
        self.assertEqual(over.right_signs(), (None, -1, 1))
        self.assertEqual(under.left_signs(), (None, -1, 1))
        self.assertEqual(under.middle_signs(), (None, 1, -1))
        self.assertEqual(under.right_signs(), (None, 1, -1))

    def test_y_pos(self):
        self.assertEqual(cup.left_strand_y_pos(1), None)
        self.assertEqual(cup.left_strand_y_pos(2), None)
        self.assertEqual(cup.middle_strand_y_pos(1), 1)
        self.assertEqual(cup.middle_strand_y_pos(2), 1)
        self.assertEqual(cup.right_strand_y_pos(1), .5)
        self.assertEqual(cup.right_strand_y_pos(2), 1.5)
        self.assertEqual(cap.left_strand_y_pos(1), .5)
        self.assertEqual(cap.left_strand_y_pos(2), 1.5)
        self.assertEqual(cap.middle_strand_y_pos(1), 1)
        self.assertEqual(cap.middle_strand_y_pos(2), 1)
        self.assertEqual(cap.right_strand_y_pos(1), None)
        self.assertEqual(cap.right_strand_y_pos(2), None)
        self.assertEqual(over.left_strand_y_pos(1), .5)
        self.assertEqual(over.left_strand_y_pos(2), 1.5)
        self.assertEqual(over.middle_strand_y_pos(1), .5)
        self.assertEqual(over.middle_strand_y_pos(2), 1.5)
        self.assertEqual(over.right_strand_y_pos(1), 1.5)
        self.assertEqual(over.right_strand_y_pos(2), .5)
        self.assertEqual(under.left_strand_y_pos(1), 1.5)
        self.assertEqual(under.left_strand_y_pos(2), .5)
        self.assertEqual(under.middle_strand_y_pos(1), .5)
        self.assertEqual(under.middle_strand_y_pos(2), 1.5)
        self.assertEqual(under.right_strand_y_pos(1), .5)
        self.assertEqual(under.right_strand_y_pos(2), 1.5)

    def test_comparison(self):
        self.assertEqual(t1 + t2, t1 + ETangle(ETangle.Type.CUP, (-1, 1, -1, 1), 3))

    def test_sum(self):
        self.assertEqual(trefoil, t1 + t2 + t3 + t4 + t5 + t6 + t7)
        self.assertEqual(trefoil, sum([t1, t2, t3, t4, t5, t6, t7]))

    def test_ascii_rendering(self):
        print('\n')
        print(unknot2)
        print(trefoil)
        print(torus_knot)

    def test_svg_rendering(self):
        TangleRenderer.svg('output/trefoil.svg', trefoil)
        TangleRenderer.svg('output/unknot2.svg', unknot2)


if __name__ == '__main__':
    unittest.main()
