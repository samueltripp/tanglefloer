from SignAlgebra.AMinus import *

# A few tests
am = AMinus([-1, 1, 1, -1])
r = am.polyring
U1 = am.polyring['U1']
U2 = am.polyring['U2']

# First multiplication of single elements
# double crosses black, should be empty product
elt1 = am.gen({4: 0, 3: 1})
elt2 = am.gen({0: 4, 1: 3})
assert elt1 * elt2 == am.zero()

# double crosses two rightward strands, should be nontrivial product
elt3 = am.gen({4: 0, 3: 1})
elt4 = am.gen({0: 0, 1: 3})
assert elt3 * elt4 == (U1 * U2) * am.gen({4: 0, 3: 3})

# just straightforward multiplication
elt5 = am.gen({4: 1, 3: 2})
elt6 = am.gen({2: 1, 1: 0})
assert elt5 * elt6 == am.gen({4: 0, 3: 1})

# testing differentials
elt = am.gen({0: 4, 1: 2, 3: 1})
assert elt.diff() == U1 * am.gen({0: 4, 1: 1, 3: 2}) + am.gen({0: 2, 1: 4, 3: 1})

# testing gradings
# TODO