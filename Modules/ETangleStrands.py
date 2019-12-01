from __future__ import annotations
from Modules.CTMinus import *
from frozendict import *
from Modules.StrandDiagram import StrandDiagram
from Tangles.Tangle import *
from SignAlgebra.AMinus import *


# represents a pair of partial bijections overlaid on an elementary tangle
class ETangleStrands:
    def __init__(self, etangle: ETangle, left_strands: Dict, right_strands: Dict):
        self.etangle = etangle
        assert validdictionaries(etangle, left_strands, right_strands)
        self.left_strands = frozendict(left_strands)
        self.left_strands_inverse = frozendict(invert_injection(left_strands))
        self.right_strands = frozendict(right_strands)
        self.right_strands_inverse = frozendict(invert_injection(right_strands))

    # the idempotent e^D_L                                                                                                                                
    def left_idempotent(self) -> AMinus.Element:
        return self.etangle.left_algebra.idempotent(self.left_idempotent_strands().keys())

    def left_idempotent_strands(self) -> Dict:
        occupied = self.left_strands.keys()
        total = set(range(len(self.etangle.left_algebra.ss)))
        return {strand: strand for strand in total - occupied}

    # the idempotent e^A_R                                                                                                                                
    def right_idempotent(self) -> AMinus.Element:
        return self.etangle.right_algebra.idempotent(list(self.right_strands.values()))

    def left_y_pos(self, black_strand: int):
        return self.left_strands_inverse[black_strand]

    def right_y_pos(self, black_strand: int):
        return self.right_strands[black_strand]

    def to_strand_diagram(self):
        orange_strands = {}
        orange_signs = {}
        for orange in range(1, len(self.etangle.signs)):
            orange_strands[orange] = (
                self.etangle.left_y_pos(orange), self.etangle.middle_y_pos(orange), self.etangle.right_y_pos(orange)
            )
            orange_signs[orange] = self.etangle.middle_signs()[orange]
        black_strands = {}
        for black in self.left_strands.values():
            black_strands[black] = (self.left_y_pos(black), black, None)
        for black in self.right_strands.keys():
            black_strands[black] = (None, black, self.right_y_pos(black))

        return StrandDiagram(orange_strands, orange_signs, black_strands)

    def left_idempotent_strand_diagram(self):
        orange_strands = {}
        orange_signs = {}
        for orange in range(1, len(self.etangle.signs)):
            if self.etangle.left_y_pos(orange):
                orange_strands[orange] = (
                    self.etangle.left_y_pos(orange), self.etangle.left_y_pos(orange), self.etangle.middle_y_pos(orange)
                )
            orange_signs[orange] = self.etangle.left_signs()[orange]
        black_strands = {}
        for black in self.etangle.left_points():
            if black in self.left_strands.keys():
                black_strands[black] = (None, black, self.left_strands[black])
            else:
                black_strands[black] = (black, black, None)
        return StrandDiagram(orange_strands, orange_signs, black_strands)

    def __repr__(self):
        return str((self.etangle, self.left_strands, self.right_strands))

    def __eq__(self, other: ETangleStrands):
        return other is not None and \
               self.etangle == other.etangle and \
               self.left_strands == other.left_strands and \
               self.right_strands == other.right_strands

    def __hash__(self):
        return hash(self.etangle) + hash(self.left_strands) + hash(self.right_strands)


def validdictionaries(et: ETangle, ls: Dict, rs: Dict):
    for key in ls.keys():
        if key not in et.left_points():
            return False
    for val in rs.values():
        if val not in et.right_points():
            return False
    return set(ls.values()).union(set(rs.keys())) == set(et.middle_points())
