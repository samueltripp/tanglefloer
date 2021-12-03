from Tangles.Tangle import *
from Tangles.TangleRenderer import *


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

def test_signs():
    assert cup.left_sign_sequence() == (None,)
    assert cup.right_sign_sequence() == (None, 1, -1)
    assert cap.left_sign_sequence() == (None, 1, -1)
    assert cap.right_sign_sequence() == (None,)
    assert over.left_sign_sequence() == (None, 1, -1)
    assert over.right_sign_sequence() == (None, -1, 1)
    assert under.left_sign_sequence() == (None, -1, 1)
    assert under.right_sign_sequence() == (None, 1, -1)

    assert cup2.left_sign_sequence() == (None, -1, -1)
    assert over2.left_sign_sequence() == (None, -1, -1, 1, -1)
    assert under2.left_sign_sequence() == (None, -1, 1, -1, -1)
    assert cap2.left_sign_sequence() == (None, -1, -1, 1, -1)
    assert cup2.right_sign_sequence() == (None, -1, -1, 1, -1)
    assert over2.right_sign_sequence() == (None, -1, 1, -1, -1)
    assert under2.right_sign_sequence() == (None, -1, -1, 1, -1)
    assert cap2.right_sign_sequence() == (None, -1, -1)

    assert [None, None] == [cup.strand_index_to_left_sign(i) for i in range(1, 3)]
    assert [1, -1] == [over.strand_index_to_left_sign(i) for i in range(1, 3)]
    assert [1, -1] == [under.strand_index_to_left_sign(i) for i in range(1, 3)]
    assert [1, -1] == [cap.strand_index_to_left_sign(i) for i in range(1, 3)]
    assert [-1, None, None, -1] == [cup2.strand_index_to_left_sign(i) for i in range(1, 5)]
    assert [-1, -1, 1, -1] == [over2.strand_index_to_left_sign(i) for i in range(1, 5)]
    assert [-1, -1, 1, -1] == [under2.strand_index_to_left_sign(i) for i in range(1, 5)]
    assert [-1, -1, 1, -1] == [cap2.strand_index_to_left_sign(i) for i in range(1, 5)]

def test_y_pos():
    assert cup.left_y_pos(1) is None
    assert cup.left_y_pos(2) is None
    assert cup.middle_y_pos(1) == 1
    assert cup.middle_y_pos(2) == 1
    assert cup.right_y_pos(1) == .5
    assert cup.right_y_pos(2) == 1.5
    assert cap.left_y_pos(1) == .5
    assert cap.left_y_pos(2) == 1.5
    assert cap.middle_y_pos(1) == 1
    assert cap.middle_y_pos(2) == 1
    assert cap.right_y_pos(1) is None
    assert cap.right_y_pos(2) is None
    assert over.left_y_pos(1) == .5
    assert over.left_y_pos(2) == 1.5
    assert over.middle_y_pos(1) == .5
    assert over.middle_y_pos(2) == 1.5
    assert over.right_y_pos(1) == 1.5
    assert over.right_y_pos(2) == .5
    assert under.left_y_pos(1) == 1.5
    assert under.left_y_pos(2) == .5
    assert under.middle_y_pos(1) == .5
    assert under.middle_y_pos(2) == 1.5
    assert under.right_y_pos(1) == .5
    assert under.right_y_pos(2) == 1.5

def test_comparison():
    assert t1 + t2 == t1 + ETangle(ETangle.Type.CUP, (-1, 1, -1, 1), 3)

def test_sum():
    assert t1 + t2 + t3 + t4 + t5 + t6 + t7 == trefoil
    assert sum([t1, t2, t3, t4, t5, t6, t7]) == trefoil

def test_ascii_rendering():
    print('\n')
    print(unknot2)
    print(trefoil)
    print(torus_knot)
    print(cup+straight+straight+cap)

def test_svg_rendering():
    TangleRenderer.svg('Test/output/trefoil.svg', trefoil)
    TangleRenderer.svg('Test/output/unknot2.svg', unknot2)
    TangleRenderer.svg('Test/output/cup.svg', cup, [{0: 2}, {0: 1}])
