from Tangles.Tangle import *
from Modules.CTMinus import *


# tests for dminus on caps
cap = ETangle(ETangle.Type.CAP, (1, 1, -1, 1, 1, 1), 3)
cap2 = ETangle(ETangle.Type.CAP, (1, 1, 1, -1, 1, 1), 3)


def fill_right(n, left):
    return {i: i for i in set(range(n)) - set(left.values())}


left1 = {1: 1, 5: 2}
sd1 = StrandDiagram(cap, left1, fill_right(7, left1))
print(sd1)
print(dminus(sd1))
