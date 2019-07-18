from Tangles.Tangle import *
from Modules.CTMinus import *


cup = ETangle(ETangle.Type.CUP, (1, -1), 1)
over = ETangle(ETangle.Type.OVER, (1, -1), 1)
under = ETangle(ETangle.Type.UNDER, (-1, 1), 1)
cap = ETangle(ETangle.Type.CAP, (1, -1), 1)
unknot = Tangle([cup, cap])
unknot2 = Tangle([cup, over, under, cap])


e1 = ETangle(ETangle.Type.OVER,(1,1,-1,-1),2)
sd = StrandDiagram(e1,{2:2,3:3,4:4},{0:1,1:0})
dplus(sd) == {StrandDiagram(e1,{2:2,3:3,4:4},{0:0,1:1}):e1.polyring['U1']}
