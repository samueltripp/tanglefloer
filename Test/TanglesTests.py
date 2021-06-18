from Tangles.Tangle import *
from Tangles.TangleRenderer import *
import unittest


# some simple examples
cup = ETangle(ETangle.Type.CUP, (1, -1), 1)
over = ETangle(ETangle.Type.OVER, (1, -1), 1)
under = ETangle(ETangle.Type.UNDER, (1, -1), 1)
cap = ETangle(ETangle.Type.CAP, (1, -1), 1)
unknot = Tangle((cup, cap))
unknot2 = Tangle((cup, over, under, cap))
straight = ETangle(ETangle.Type.STRAIGHT, (1, -1))
straight2 = ETangle(ETangle.Type.STRAIGHT, (1, 1))

# some larger examples
cup2 = ETangle(ETangle.Type.CUP, (-1, -1, 1, -1), 2)
over2 = ETangle(ETangle.Type.OVER, (-1, -1, 1, -1), 2)
under2 = ETangle(ETangle.Type.UNDER, (-1, -1, 1, -1), 2)
cap2 = ETangle(ETangle.Type.CAP, (-1, -1, 1, -1), 2)

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

    def trefoil_test(self):
        tref_da = reduced_type_da(trefoil)
        tref_cc = tref_da.to_chain_complex()
        print(tref_cc.d_squared_is_zero())
#        tref_da.to_agraph().draw('Test/output/tref_da.svg')
#        tref_cc.to_agraph().draw('Test/output/tref_cc.svg')
        tref_cc.write_m2_def('Test/output/tref_cc.m2')
        
    def test_signs(self):
        self.assertEqual(cup.left_sign_sequence(), (None,))
        self.assertEqual(cup.right_sign_sequence(), (None, 1, -1))
        self.assertEqual(cap.left_sign_sequence(), (None, 1, -1))
        self.assertEqual(cap.right_sign_sequence(), (None,))
        self.assertEqual(over.left_sign_sequence(), (None, 1, -1))
        self.assertEqual(over.right_sign_sequence(), (None, -1, 1))
        self.assertEqual(under.left_sign_sequence(), (None, -1, 1))
        self.assertEqual(under.right_sign_sequence(), (None, 1, -1))

        self.assertEqual(cup2.left_sign_sequence(), (None, -1, -1))
        self.assertEqual(over2.left_sign_sequence(), (None, -1, -1, 1, -1))
        self.assertEqual(under2.left_sign_sequence(), (None, -1, 1, -1, -1))
        self.assertEqual(cap2.left_sign_sequence(), (None, -1, -1, 1, -1))
        self.assertEqual(cup2.right_sign_sequence(), (None, -1, -1, 1, -1))
        self.assertEqual(over2.right_sign_sequence(), (None, -1, 1, -1, -1))
        self.assertEqual(under2.right_sign_sequence(), (None, -1, -1, 1, -1))
        self.assertEqual(cap2.right_sign_sequence(), (None, -1, -1))

        self.assertEqual([None, None], [cup.strand_index_to_left_sign(i) for i in range(1, 3)])
        self.assertEqual([1, -1], [over.strand_index_to_left_sign(i) for i in range(1, 3)])
        self.assertEqual([1, -1], [under.strand_index_to_left_sign(i) for i in range(1, 3)])
        self.assertEqual([1, -1], [cap.strand_index_to_left_sign(i) for i in range(1, 3)])
        self.assertEqual([-1, None, None, -1], [cup2.strand_index_to_left_sign(i) for i in range(1, 5)])
        self.assertEqual([-1, -1, 1, -1], [over2.strand_index_to_left_sign(i) for i in range(1, 5)])
        self.assertEqual([-1, -1, 1, -1], [under2.strand_index_to_left_sign(i) for i in range(1, 5)])
        self.assertEqual([-1, -1, 1, -1], [cap2.strand_index_to_left_sign(i) for i in range(1, 5)])

    def test_y_pos(self):
        self.assertEqual(cup.left_y_pos(1), None)
        self.assertEqual(cup.left_y_pos(2), None)
        self.assertEqual(cup.middle_y_pos(1), 1)
        self.assertEqual(cup.middle_y_pos(2), 1)
        self.assertEqual(cup.right_y_pos(1), .5)
        self.assertEqual(cup.right_y_pos(2), 1.5)
        self.assertEqual(cap.left_y_pos(1), .5)
        self.assertEqual(cap.left_y_pos(2), 1.5)
        self.assertEqual(cap.middle_y_pos(1), 1)
        self.assertEqual(cap.middle_y_pos(2), 1)
        self.assertEqual(cap.right_y_pos(1), None)
        self.assertEqual(cap.right_y_pos(2), None)
        self.assertEqual(over.left_y_pos(1), .5)
        self.assertEqual(over.left_y_pos(2), 1.5)
        self.assertEqual(over.middle_y_pos(1), .5)
        self.assertEqual(over.middle_y_pos(2), 1.5)
        self.assertEqual(over.right_y_pos(1), 1.5)
        self.assertEqual(over.right_y_pos(2), .5)
        self.assertEqual(under.left_y_pos(1), 1.5)
        self.assertEqual(under.left_y_pos(2), .5)
        self.assertEqual(under.middle_y_pos(1), .5)
        self.assertEqual(under.middle_y_pos(2), 1.5)
        self.assertEqual(under.right_y_pos(1), .5)
        self.assertEqual(under.right_y_pos(2), 1.5)

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
        print(cup+straight+straight+cap)

    def test_svg_rendering(self):
        TangleRenderer.svg('Test/output/trefoil.svg', trefoil)
        TangleRenderer.svg('Test/output/unknot2.svg', unknot2)
        TangleRenderer.svg('Test/output/cup.svg', cup, [{0: 2}, {0: 1}])


if __name__ == '__main__':
    unittest.main()
