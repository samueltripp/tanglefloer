from SignAlgebra.AMinus import *

# A few tests

# First multiplication of single elements
# double crosses black, should be empty product
am = AMinus([-1, 1, 1, -1])
gena = AMinusGen(am, {4: 0, 3: 1})
genb = AMinusGen(am, {0: 4, 1: 3})
assert gena * genb == am.zero()
assert gena * genb == gena.to_element() * genb.to_element()

# double crosses two rightward strands, should be nontrivial product
gena = AMinusGen(am, {4: 0, 3: 1})
genb = AMinusGen(am, {0: 0, 1: 3})
assert gena * genb == AMinusElement(am, {AMinusGen(am, {4: 0, 3: 3}): am.polyring['U1'] * am.polyring['U2']})
assert gena * genb == AMinusElement(am, {gena: 1}) * AMinusElement(am, {genb: 1})

# just straightforward multiplication
gena = AMinusGen(am, {4: 1, 3: 2})
genb = AMinusGen(am, {2: 1, 1: 0})
assert gena * genb == AMinusGen(am, {4: 0, 3: 1}).to_element()

# testing differentials
dgen = AMinusGen(am, {0: 4, 1: 2, 3: 1})
assert am.diff(dgen) == \
       AMinusElement(am, {AMinusGen(am, {0: 4, 1: 1, 3: 2}): am.polyring['U1'], AMinusGen(am, {0: 2, 1: 4, 3: 1}): 1})
