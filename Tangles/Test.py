from Tangles import *

cup = ETangle(ETangle.Type.CUP, (1,-1), 1)
over = ETangle(ETangle.Type.OVER, (1,-1), 1)
under = ETangle(ETangle.Type.UNDER, (-1,1), 1)
cap = ETangle(ETangle.Type.CAP, (1,-1), 1)
s = Tangle([cup, cap])
t = Tangle([cup, over, under, cap])
print(t.draw_ascii())

s=Tangle([ETangle(ETangle.Type.CUP, (1,-1), 1),ETangle(ETangle.Type.CUP, (1,-1,1,-1),1),ETangle(ETangle.Type.CUP, (1,-1,1,-1,1,-1),1),ETangle(ETangle.Type.CAP, (1,-1,1,-1,1,-1),1),ETangle(ETangle.Type.CAP, (1,-1,1,-1),1),ETangle(ETangle.Type.CAP, (1,-1), 1)])
print(s.draw_ascii())

t1 = ETangle(ETangle.Type.CUP, (-1,1), 1)
t2 = ETangle(ETangle.Type.CUP, (-1,1,-1,1), 3)
t3 = ETangle(ETangle.Type.OVER, (-1,1,-1,1), 2)
t4 = ETangle(ETangle.Type.UNDER, (-1,-1,1,1), 1)
t5 = ETangle(ETangle.Type.OVER, (-1,-1,1,1), 2)
t6 = ETangle(ETangle.Type.CAP, (-1,1,-1,1), 1)
t7 = ETangle(ETangle.Type.CAP, (-1,1), 1)
t = Tangle([t1,t2,t3,t4,t5,t6,t7])
print(t.draw_ascii())