from __future__ import annotations
from Tangles.Tangle import *
from Modules.Bimodule import *
from Tangles.Functions import *


# represents a pair of partial bijections overlaid on an elementary tangle
class StrandDiagram:
    def __init__(self, etangle: ETangle, left_strands: Dict, right_strands: Dict):
        self.etangle = etangle
        self.left_strands = left_strands
        self.right_strands = right_strands

    def __repr__(self):
        return str((self.etangle, self.left_strands, self.right_strands))


def type_da(etangle: ETangle) -> Bimodule:
    left_algebra = AMinus(etangle.left_signs())
    right_algebra = AMinus(etangle.right_signs())
    gens = [StrandDiagram(etangle, left_strands, right_strands) for left_strands, right_strands in
            enumerate_gens([etangle.left_points(), etangle.middle_points(), etangle.right_points()])]
    maps = None  # TODO

    return Bimodule(left_algebra, right_algebra, gens, maps)


def reverse_injection(d):
    return {v: k for k, v in d.items()}


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
