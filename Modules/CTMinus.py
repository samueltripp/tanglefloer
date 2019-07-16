from __future__ import annotations
from typing import List
from Tangles.Tangle import *
from Modules.Bimodule import *
from Tangles.Functions import *


# represents a pair of partial bijections overlaid on an elementary tangle
class StrandDiagram:
    def __init__(self, etangle: ETangle, left_strands: Dict, right_strands: Dict):
        self.etangle = etangle
        self.left_strands = left_strands
        self.right_strands = right_strands

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


def type_da(etangle: ETangle) -> Bimodule:
    strand_diagrams = [StrandDiagram(etangle, left_strands, right_strands)
                       for left_strands, right_strands in
                       enumerate_gens([etangle.left_points(), etangle.middle_points(), etangle.right_points()])]
    maps = sum((delta1_1(x) for x in strand_diagrams), []) + \
        sum((delta1_2(x, a) for x in strand_diagrams
             for a in etangle.right_algebra.left_gens(list(x.left_strands.keys()))), [])

    return Bimodule.from_strand_diagrams(etangle.left_algebra, etangle.right_algebra, strand_diagrams, maps)


def delta1_1(x: StrandDiagram) -> List[Bimodule.Edge]:
    return []  # TODO


def delta1_2(x: StrandDiagram, a: AMinusElement) -> List[Bimodule.Edge]:
    return []  # TODO


# points - a list of sets of points
def enumerate_gens(points):
    sequences = []
    if len(points) < 2:
        return sequences
    elif len(points) == 2:
        return [[pb] for pb in partial_bijections(points[0], points[1])]
    else:
        for pb in partial_bijections(points[0], points[1]):
            coker = set(points[1]).difference(pb.values())
            sequences.extend(
                [[pb] + sequence for sequence in enumerate_gens_helper([list(coker)] + points[2:])])
    return sequences


def enumerate_gens_helper(points):
    sequences = []
    if len(points) < 2:
        return sequences
    elif len(points) == 2:
        return [[inj] for inj in injections(points[0], points[1])]
    else:
        for inj in injections(points[0], points[1]):
            coker = set(points[1]).difference(inj.values())
            sequences.extend(
                [[inj] + sequence for sequence in enumerate_gens_helper([list(coker)] + points[2:])])
    return sequences
