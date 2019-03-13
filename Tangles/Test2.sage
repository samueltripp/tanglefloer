#load("Tangle.sage")
#load("EnhancedTangle.sage")

# some example tangles
cup = ETangle(ETangle.Type.CUP, (1,-1), 1)
over = ETangle(ETangle.Type.OVER, (1,1), 1)
under = ETangle(ETangle.Type.UNDER, (1,-1), 1)
cap = ETangle(ETangle.Type.CAP, (1,-1), 1)

Ecup = EnhancedETangle(cup, {0:2}, {0:1})
Eover = EnhancedETangle(over, {1:1}, {0:1,2:0})
Eunder = EnhancedETangle(under, {1:1}, {0:1, 2:0})


A = Eunder.Dplus()
print A
#for a in A.keys():
#    print a.bijplus
#    print A[a]
print Eover.Dplus()
