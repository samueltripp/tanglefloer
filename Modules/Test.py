from Tangles.Tangle import *
from Modules.CTMinus import *


# tests for dminus on caps
cap = ETangle(ETangle.Type.CAP, (1, 1, -1, 1, 1, 1), 3)
cup = ETangle(ETangle.Type.CUP, (1, 1, -1, 1, 1, 1), 3)
over = ETangle(ETangle.Type.OVER, (1, 1, -1, 1, 1, 1), 3)
under = ETangle(ETangle.Type.UNDER, (1, 1, -1, 1, 1, 1), 3)


lefts = [{1: 1, 5: 2}, {3: 2, 5: 3}, {4: 3, 6: 4}, {1: 1, 5: 3}, {0: 0, 1: 4}, {1: 1, 3: 3}]


def fill_right(etangle, left_strands):
    out = {}
    right_points = set(etangle.right_points())
    for p in set(etangle.middle_points()) - set(left_strands.values()):
        out[p] = right_points.pop()
    return out


print('testing cap')
for left in lefts:
    sd = StrandDiagram(cap, left, fill_right(cap, left))
    print('strand diagram: ' + str(sd))
    print('dminus: ' + str(dminus(sd)))
print('\n')
print('testing cup')
for left in lefts:
    sd = StrandDiagram(cup, left, fill_right(cup, left))
    print('strand diagram: ' + str(sd))
    print('dminus: ' + str(dminus(sd)))
print('\n')
print('testing over')
for left in lefts:
    sd = StrandDiagram(over, left, fill_right(over, left))
    print('strand diagram: ' + str(sd))
    print('dminus: ' + str(dminus(sd)))
print('\n')
print('testing under')
for left in lefts:
    sd = StrandDiagram(under, left, fill_right(under, left))
    print('strand diagram: ' + str(sd))
    print('dminus: ' + str(dminus(sd)))
print('\n')


