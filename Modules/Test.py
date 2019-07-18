from Tangles.Tangle import *
from Modules.CTMinus import *


# tests for dminus on caps
cap1 = ETangle(ETangle.Type.CAP, (1, 1, -1, 1, 1, 1), 3)
cap2 = ETangle(ETangle.Type.CAP, (1, 1, 1, -1, 1, 1), 3)


def fill_right(n, left_strands):
    return {i: i for i in set(range(n)) - set(left_strands.values())}


lefts = [{1: 1, 5: 2}, {3: 2, 5: 4}, {4: 4, 6: 5}, {1: 1, 5: 4}, {0: 0, 1: 5}, {1: 1, 3: 4}]

print('cap1')
for left in lefts:
    sd = StrandDiagram(cap1, left, fill_right(7, left))
    print('strand diagram: ' + str(sd))
    print('dminus: ' + str(dminus(sd)))

print('cap2')
for left in lefts:
    sd = StrandDiagram(cap2, left, fill_right(7, left))
    print('strand diagram: ' + str(sd))
    print('dminus: ' + str(dminus(sd)))
