from __future__ import annotations
from Modules.CTMinus import *
from frozendict import *
from Tangles.Tangle import *
from SignAlgebra.AMinus import *


# represents a pair of partial bijections overlaid on an elementary tangle
class StrandDiagram:
    def __init__(self, etangle: ETangle, left_strands: Dict, right_strands: Dict):
        self.etangle = etangle
        self.left_strands = frozendict(left_strands)
        self.right_strands = frozendict(right_strands)

    # the idempotent e^D_L                                                                                                                                
    def left_idempotent(self):
        occupied = self.left_strands.keys()
        total = set(range(len(self.etangle.left_algebra.ss)))
        return self.etangle.left_algebra.idempotent(list(total - occupied))

    # the idempotent e^A_R                                                                                                                                
    def right_idempotent(self) -> AMinusElement:
        return self.etangle.right_algebra.idempotent(list(self.right_strands.values()))

    def __repr__(self):
        return str((self.etangle, self.left_strands, self.right_strands))

    def __eq__(self, other: StrandDiagram):
        return self.etangle == other.etangle and \
                self.left_strands == other.left_strands and \
                self.right_strands == other.right_strands

    def __hash__(self):
        return hash(self.etangle) + hash(self.left_strands) + hash(self.right_strands)
