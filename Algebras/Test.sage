load('Algebras/Generator.sage')
load('Algebras/SignSequence.sage')

# These are the examples in Figure 2 of Intro To Tangle Floer
a = SignSequence([-1,1,1,-1])
b = Generator(a,[4,2,-1,1,-1])
c = b.differential()
print c.dict

d = Generator(a,[-1,-1,-1,2,1])
e = Generator(a,[-1,0,1,-1,-1])
print (d*e).dict
