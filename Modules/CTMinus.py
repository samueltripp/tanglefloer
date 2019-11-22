from __future__ import annotations
from typing import List
from Tangles.Tangle import *
from Modules.Bimodule import *
from Modules.Bimodule import Bimodule
from Tangles.Functions import *
from Modules.StrandDiagram import *


def type_da(etangle: ETangle) -> TypeDA:
    strand_diagrams = [StrandDiagram(etangle, left_strands, right_strands)
                       for left_strands, right_strands in
                       enumerate_gens([etangle.left_points(), etangle.middle_points(), etangle.right_points()])]
    maps = sum((delta1_1(x) for x in strand_diagrams), []) + \
           [delta1_2(x, a) for x in strand_diagrams
            for a in etangle.right_algebra.left_gens(list(x.left_strands.keys()))]

    return TypeDA.from_strand_diagrams(etangle.left_algebra, etangle.right_algebra, strand_diagrams, maps)


def delta1_1(x: StrandDiagram) -> List[Bimodule.Edge]:
    out = []
    out += [Bimodule.Edge(x, y, c, (x.left_idempotent(),), tuple()) for y, c in dplus(x).d.items()]
    out += [Bimodule.Edge(x, y, c, (x.left_idempotent(),), tuple()) for y, c in dminus(x).d.items()]
    out += [Bimodule.Edge(x, y, c, (x.left_idempotent(),), tuple()) for y, c in dmixed(x).d.items()]
    out += [delta_ell(x)]
    return out


def delta1_2(x: StrandDiagram, a: AMinusElement) -> Bimodule.Edge:
    e = m2(x, a)
    return Bimodule.Edge(x, e.target_diagram, e.c, (x.left_idempotent(),), (a,))


def m2(x: StrandDiagram, a: AMinusElement) -> Bimodule.Edge:
    pass  # TODO


def delta_ell(x: StrandDiagram) -> Bimodule.Edge:
    pass  # TODO


def dplus(sd: StrandDiagram) -> Bimodule.Element:
    out = Bimodule.Element()
    for p1 in sd.right_strands.keys():
        for p2 in sd.right_strands.keys():
            # if the black strands cross
            if p1 > p2 and sd.right_strands[p1] < sd.right_strands[p2]:
                # smooth the crossing and add it to the output
                out += smooth_right_crossing(sd, p1, p2)
    return out


def smooth_right_crossing(sd: StrandDiagram, p1: int, p2: int) -> Bimodule.Element:
    c = sd.etangle.polyring.one()
    q1 = sd.right_strands[p1]
    q2 = sd.right_strands[p2]
    s = sd.etangle.position
    sd_out = StrandDiagram(sd.etangle, sd.left_strands, swap_values(sd.right_strands, p1, p2))

    # first, check for black-black double-crossings
    for p3 in range(p2+1, p1):
        if p3 in sd.right_strands:
            q3 = sd.right_strands[p3]
            if q1 < q3 < q2:
                return Bimodule.Element()

    # next, check for black-orange double-crossings
    orange_strands_double_crossed = []

    # for each orange strand, we want to check if it double-crosses either of the black strands
    for orange_strand in range(1, len(sd.etangle.signs)):
        right_y_pos = sd.etangle.right_strand_y_pos(orange_strand)
        if right_y_pos is None:
            continue
        middle_y_pos = sd.etangle.middle_strand_y_pos(orange_strand)
        times_crossed_p1 = 0
        times_crossed_p2 = 0

        # count how many times this orange strand crosses the p1 -> q2 black strand
        if (p1 < middle_y_pos) ^ (min(p1, q2) < middle_y_pos):
            times_crossed_p1 += 1
        if (min(p1, q2) < middle_y_pos) ^ (q2 < middle_y_pos):
            times_crossed_p1 += 1
        if (q2 < middle_y_pos) ^ (q2 < right_y_pos):
            times_crossed_p1 += 1

        # count how many times this orange strand crosses the p2 -> q1 black strand
        if (p2 < middle_y_pos) ^ (min(p1, q2) < middle_y_pos):
            times_crossed_p2 += 1
        if (min(p1, q2) < middle_y_pos) ^ (q1 < middle_y_pos):
            times_crossed_p2 += 1
        if (q1 < middle_y_pos) ^ (q1 < right_y_pos):
            times_crossed_p2 += 1

        # if either black strand is double-crossed, add this orange strand to the list
        if times_crossed_p1 > 1 or times_crossed_p2 > 1:
            orange_strands_double_crossed += [orange_strand]

    # turn the list of orange strands into a coefficient
    for orange_strand in orange_strands_double_crossed:
        if sd.etangle.middle_signs()[orange_strand] == 1:
            c *= sd.etangle.polyring['U' + str(orange_strand)]
        else:
            return Bimodule.Element()

    return Bimodule.Element({sd_out: c})


def dminus(sd: StrandDiagram) -> Bimodule.Element:
    out = Bimodule.Element()
    for p1 in sd.left_strands.keys():
        for p2 in sd.left_strands.keys():
            # if the black strands don't cross
            if p1 > p2 and sd.left_strands[p1] > sd.left_strands[p2]:
                # introduce a crossing and add it to the output
                out += introduce_left_crossing(sd, p1, p2)
    return out


def introduce_left_crossing(sd: StrandDiagram, p1: int, p2: int) -> Bimodule.Element:
    c = sd.etangle.polyring.one()
    q1 = sd.left_strands[p1]
    q2 = sd.left_strands[p2]
    s = sd.etangle.position
    sd_out = StrandDiagram(sd.etangle, swap_values(sd.left_strands, p1, p2), sd.right_strands)

    # first, check for black-black double-crossings
    for p3 in range(p2+1, p1):
        if p3 in sd.left_strands:
            q3 = sd.left_strands[p3]
            if q2 < q3 < q1:
                return Bimodule.Element()

    # next, check for black-orange double-crossings
    orange_strands_double_crossed = []

    # for each orange strand, we want to check if it double-crosses either of the black strands
    for orange_strand in range(1, len(sd.etangle.signs)):
        left_y_pos = sd.etangle.left_strand_y_pos(orange_strand)
        if left_y_pos is None:
            continue
        middle_y_pos = sd.etangle.middle_strand_y_pos(orange_strand)
        times_crossed_p1 = 0
        times_crossed_p2 = 0

        # count how many times this orange strand crosses the p1 -> q2 black strand
        if (p1 < left_y_pos) ^ (min(p2, q1) < left_y_pos):
            times_crossed_p1 += 1
        if (min(p2, q1) < left_y_pos) ^ (q1 < left_y_pos):
            times_crossed_p1 += 1
        if (q1 < left_y_pos) ^ (q1 < middle_y_pos):
            times_crossed_p1 += 1

        # count how many times this orange strand crosses the p2 -> q1 black strand
        if (p2 < left_y_pos) ^ (min(p2, q1) < left_y_pos):
            times_crossed_p2 += 1
        if (min(p2, q1) < left_y_pos) ^ (q2 < left_y_pos):
            times_crossed_p2 += 1
        if (q2 < left_y_pos) ^ (q2 < middle_y_pos):
            times_crossed_p2 += 1

        # if either black strand is double-crossed, add this orange strand to the list
        if times_crossed_p1 > 1 or times_crossed_p2 > 1:
            orange_strands_double_crossed += [orange_strand]

    # turn the list of orange strands into a coefficient
    for orange_strand in orange_strands_double_crossed:
        if sd.etangle.middle_signs()[orange_strand] == -1:
            c *= sd.etangle.polyring['U' + str(orange_strand)]
        else:
            return Bimodule.Element()

    return Bimodule.Element({sd_out: c})


# swap the values associated to the given keys
def swap_values(d: Dict, k1, k2) -> Dict:
    d_out = dict(d)
    v1 = d[k1]
    v2 = d[k2]
    d_out[k1] = v2
    d_out[k2] = v1
    return d_out


def dmixed(sd: StrandDiagram) -> Bimodule.Element:
    return Bimodule.Element()  # TODO


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
