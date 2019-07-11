load('Algebras/Generator.sage')
load('Algebras/SignSequence.sage')
load('Algebras/AlgElement.sage')
load('Algebras/AMinus.sage')

# These are the examples in Figure 2 of Intro To Tangle Floer
a = SignSequence([-1,1,1,-1])
b = Generator(a,[4,2,-1,1,-1])
c = b.differential()
print c.dict

d = Generator(a,[-1,-1,-1,2,1])
e = Generator(a,[-1,0,1,-1,-1])
print (d*e).dict

am = AMinus([-1,1,1,-1])
g1 = Generator(SignSequence([-1,1,1,-1]),[-1,-1,-1,2,1])
g2 = Generator(SignSequence([-1,1,1,-1]),[-1,0,1,-1,-1])
prod =  am.multiply_in_algebra(g1,g2)
print prod
