from SignAlgebra.AMinus import AMinus
from SignAlgebra.TensorAlgebra import TensorAlgebra


def test_mul():
    am = AMinus([-1, 1, 1, -1])
    r = am.ring
    u1 = am.ring['U1']
    u2 = am.ring['U2']

    # First multiplication of single elements
    # double crosses black, should be empty product
    elt1 = am.generator({4: 0, 3: 1})
    elt2 = am.generator({0: 4, 1: 3})
    assert elt1 * elt2 == am.zero()

    # double crosses two rightward strands, should be nontrivial product
    elt3 = am.generator({4: 0, 3: 1})
    elt4 = am.generator({0: 0, 1: 3})
    assert elt3 * elt4 == (u1 * u2) * am.generator({4: 0, 3: 3})

    # just straightforward multiplication
    elt5 = am.generator({4: 1, 3: 2})
    elt6 = am.generator({2: 1, 1: 0})
    assert elt5 * elt6 == am.generator({4: 0, 3: 1})


def test_diff():
    am = AMinus([-1, 1, 1, -1])
    r = am.ring
    u1 = am.ring['U1']
    u2 = am.ring['U2']

    elt = am.generator({0: 4, 1: 2, 3: 1})
    assert elt.diff() == u1 * am.generator({0: 4, 1: 1, 3: 2}) + am.generator({0: 2, 1: 4, 3: 1})


def test_gradings():
    pass  # TODO


def test_tensor_algebra():
    am = AMinus([1])
    ta = TensorAlgebra(am)
    r = am.ring
    u1 = am.ring['U1']
    a = am.generator({0: 0, 1: 1})
    b = am.generator({0: 1, 1: 0})

    z = ta.zero()
    o = ta.one()
    assert z + (a ** o) ** b == a ** (z + o ** b)
