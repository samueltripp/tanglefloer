from Tangles.Tangle import *
from Modules.CTMinus import *


# tests for dminus on caps
cap = ETangle(ETangle.Type.CAP, (1, 1, -1, 1, 1, 1), 3)
cup = ETangle(ETangle.Type.CUP, (1, 1, -1, 1, 1, 1), 3)
over = ETangle(ETangle.Type.OVER, (1, 1, -1, 1, 1, 1), 3)
under = ETangle(ETangle.Type.UNDER, (1, 1, -1, 1, 1, 1), 3)

left = [{1: 1, 5: 2}, {3: 2, 5: 3}, {4: 3, 6: 4}, {1: 1, 5: 3}, {0: 0, 1: 4}, {1: 1, 3: 3}]
cup_left = [{1: 1, 4: 2}, {3: 2, 4: 3}, {2: 1, 4: 2}, {1: 1, 4: 3}, {0: 0, 1: 4}, {1: 1, 3: 3}]

tests = {cap: left, cup: cup_left, over: left, under: left}


def fill_right(etangle, left_strands):
    out = {}
    right_points = set(etangle.right_points())
    for p in set(etangle.middle_points()) - set(left_strands.values()):
        out[p] = right_points.pop()
    return out


for etangle, strands_list in tests.items():
    print('testing ' + str(etangle.etype))
    for strands in strands_list:
        sd = StrandDiagram(etangle, strands, fill_right(etangle, strands))
        print('strand diagram: ' + str(sd))
        print('dminus: ' + str(dminus(sd)))
    print('\n')
