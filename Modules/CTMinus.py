from __future__ import annotations
from Tangles.Tangle import *
from Modules.Bimodule import *
from Tangles.Functions import *


def type_da(etangle: ETangle) -> Bimodule:
    left_algebra = AMinus(etangle.left_signs())
    right_algebra = AMinus(etangle.right_signs())
    gens = enumerate_gens([etangle.left_points(), etangle.middle_points(), etangle.right_points()])
    maps = None  # TODO

    return Bimodule(left_algebra, right_algebra, gens, maps)


# points - a list of sets of points
def enumerate_gens(points):
    def reverse_injection(d):
        return {v: k for k, v in d.items()}

    sequences = []
    if len(points) < 2:
        return sequences
    elif len(points) == 2:
        return [[pb] for pb in partial_bijections(points[0], points[1])]
    else:
        mid = len(points) // 2
        for pb in partial_bijections(points[mid], points[mid + 1]):
            r, s = list(reversed(points[:mid + 1])), points[mid + 1:]

            r[0] = list(set(r[0]).difference(pb.keys()))
            s[0] = list(set(s[0]).difference(pb.values()))

            t = list(enumerate_gens_helper(s))
            sequences.extend([[reverse_injection(inj) for inj in reversed(s1)]
                              + [pb] + s2 for s1 in enumerate_gens_helper(r) for s2 in t])
        return sequences


def enumerate_gens_helper(points):
    if len(points) == 2:
        for inj in injections(points[0], points[1]):
            yield [inj]
    elif len(points) > 2:
        for inj in injections(points[0], points[1]):
            coker = set(points[1]).difference(inj.values())
            for sequence in enumerate_gens_helper([list(coker)] + points[2:]):
                yield [inj] + sequence
