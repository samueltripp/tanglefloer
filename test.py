load("AMinus.py")
load("AMinusGen.py")

# A few tests

# First multiplication of single elements
# double crosses black, should be empty product
am = AMinus([-1,1,1,-1])
gena = AMinusGen({4:0,3:1})
genb = AMinusGen({0:4,1:3})
assert am.multiply(gena,genb)=={}
assert cmp(am.multiply(gena,genb),am.multiply({gena:1},{genb:1})) == 0

# double crosses two rightward strands, should be nontrivial product
gena = AMinusGen({4:0,3:1})
genb = AMinusGen({0:0,1:3})
assert cmp(am.multiply(gena,genb),{AMinusGen({4:0,3:3}):am.polyring('U1')*am.polyring('U2')}) == 0
assert cmp(am.multiply(gena,genb),am.multiply({gena:1},{genb:1}))==0

# just straightforward multiplication
gena = AMinusGen({4:1,3:2})
genb = AMinusGen({2:1,1:0})
assert cmp(am.multiply(gena,genb),{AMinusGen({4:0,3:1}):am.polyring(1)})==0

# testing differentials
dgen = AMinusGen({0:4,1:2,3:1})
assert cmp(am.diff(dgen),{AMinusGen({0:4,1:1,3:2}):am.polyring('U1'),AMinusGen({0:2,1:4,3:1}):am.polyring(1)})==0
