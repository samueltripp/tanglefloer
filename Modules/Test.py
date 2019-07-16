from Tangles.Tangle import *
from Modules.CTMinus import *


cup = ETangle(ETangle.Type.CUP, (1, -1), 1)
over = ETangle(ETangle.Type.OVER, (1, -1), 1)
under = ETangle(ETangle.Type.UNDER, (-1, 1), 1)
cap = ETangle(ETangle.Type.CAP, (1, -1), 1)
unknot = Tangle([cup, cap])
unknot2 = Tangle([cup, over, under, cap])

print(enumerate_gens([over.left_points(), over.middle_points(), over.right_points()]))
