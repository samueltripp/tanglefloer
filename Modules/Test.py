from Tangles.Tangle import *
from Modules.CTMinus import *


# tests for dminus on caps
cap1 = ETangle(ETangle.Type.CAP, (1, 1, -1, 1, 1, 1), 3)
cap2 = ETangle(ETangle.Type.CAP, (1, 1, 1, -1, 1, 1), 3)


def fill_right(etangle, left_strands):
    out = {}
    right_points = set(etangle.right_points())
    for p in set(etangle.middle_points()) - set(left_strands.values()):
        out[p] = right_points.pop()
    return out


lefts = [{1: 1, 5: 2}, {3: 2, 5: 3}, {4: 3, 6: 4}, {1: 1, 5: 3}, {0: 0, 1: 4}, {1: 1, 3: 3}]

print('cap1')
for left in lefts:
    sd = StrandDiagram(cap1, left, fill_right(cap1, left))
    print('strand diagram: ' + str(sd))
    print('dminus: ' + str(dminus(sd)))

print('cap2')
for left in lefts:
    sd = StrandDiagram(cap2, left, fill_right(cap1, left))
    print('strand diagram: ' + str(sd))
    print('dminus: ' + str(dminus(sd)))
