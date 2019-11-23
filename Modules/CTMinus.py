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


def dplus(sd: StrandDiagram) -> Bimodule.Element:
    out = Bimodule.Element()
    for b1 in sd.right_strands.keys():
        for b2 in sd.right_strands.keys():
            # if the black strands cross
            if b1 > b2 and sd.right_y_pos(b1) < sd.right_y_pos(b2):
                # smooth the crossing and add it to the output
                out += smooth_right_crossing(sd, b1, b2)
    return out


def dminus(sd: StrandDiagram) -> Bimodule.Element:
    out = Bimodule.Element()
    for b1 in sd.left_strands.values():
        for b2 in sd.left_strands.values():
            # if the black strands don't cross
            if b1 > b2 and sd.left_y_pos(b1) > sd.left_y_pos(b2):
                # introduce a crossing and add it to the output
                out += introduce_left_crossing(sd, b1, b2)
    return out


def dmixed(sd: StrandDiagram) -> Bimodule.Element:
    out = Bimodule.Element()

    for b1 in sd.right_strands.keys():
        for b2 in sd.right_strands.keys():
            # if two black strands don't cross on the right
            if b1 > b2 and sd.right_y_pos(b1) > sd.right_y_pos(b2):
                out += dmixed_case1(sd, b1, b2)

    for b1 in sd.left_strands.values():
        for b2 in sd.left_strands.values():
            # if two black strands cross on the left
            if b1 > b2 and sd.left_y_pos(b1) < sd.left_y_pos(b2):
                out += dmixed_case2(sd, b1, b2)

    for b1 in sd.left_strands.values():
        for b2 in sd.right_strands.keys():
            if b1 > b2:
                out += dmixed_case3(sd, b1, b2)

    for b1 in sd.right_strands.keys():
        for b2 in sd.left_strands.values():
            if b1 > b2:
                out += dmixed_case4(sd, b1, b2)

    return out


def delta_ell(x: StrandDiagram) -> Bimodule.Edge:
    pass  # TODO


def m2(x: StrandDiagram, a: AMinusElement) -> Bimodule.Edge:
    pass  # TODO


def smooth_right_crossing(sd: StrandDiagram, b1: int, b2: int) -> Bimodule.Element:
    c = sd.etangle.polyring.one()
    a1 = sd.right_y_pos(b1)
    a2 = sd.right_y_pos(b2)
    sd_out = StrandDiagram(sd.etangle, sd.left_strands, swap_values(sd.right_strands, b1, b2))

    # first, check for black-black double-crossings
    for b3 in range(b2 + 1, b1):
        if b3 in sd.right_strands:
            a3 = sd.right_strands[b3]
            if a1 < a3 < a2:
                return Bimodule.Element()

    # next, check for black-orange double-crossings
    orange_strands_double_crossed = []

    # for each orange strand, we want to check if it double-crosses either of the black strands
    for orange in range(1, len(sd.etangle.signs)):
        right_y_pos = sd.etangle.right_y_pos(orange)
        if right_y_pos is None:
            continue
        middle_y_pos = sd.etangle.middle_y_pos(orange)
        times_crossed_p1 = 0
        times_crossed_p2 = 0

        # count how many times this orange strand crosses the p1 -> q2 black strand
        if (b1 < middle_y_pos) ^ (min(b1, a2) < middle_y_pos):
            times_crossed_p1 += 1
        if (min(b1, a2) < middle_y_pos) ^ (a2 < middle_y_pos):
            times_crossed_p1 += 1
        if (a2 < middle_y_pos) ^ (a2 < right_y_pos):
            times_crossed_p1 += 1

        # count how many times this orange strand crosses the p2 -> q1 black strand
        if (b2 < middle_y_pos) ^ (min(b1, a2) < middle_y_pos):
            times_crossed_p2 += 1
        if (min(b1, a2) < middle_y_pos) ^ (a1 < middle_y_pos):
            times_crossed_p2 += 1
        if (a1 < middle_y_pos) ^ (a1 < right_y_pos):
            times_crossed_p2 += 1

        # if either black strand is double-crossed, add this orange strand to the list
        if times_crossed_p1 > 1 or times_crossed_p2 > 1:
            orange_strands_double_crossed += [orange]

    # turn the list of orange strands into a coefficient
    for orange in orange_strands_double_crossed:
        if sd.etangle.middle_signs()[orange] == 1:
            c *= sd.etangle.strand_index_to_variable(orange)
        else:
            return Bimodule.Element()

    return Bimodule.Element({sd_out: c})


def introduce_left_crossing(sd: StrandDiagram, b1: int, b2: int) -> Bimodule.Element:
    c = sd.etangle.polyring.one()
    a1 = sd.left_y_pos(b1)
    a2 = sd.left_y_pos(b2)
    sd_out = StrandDiagram(sd.etangle, swap_values(sd.left_strands, a1, a2), sd.right_strands)

    # first, check for black-black double-crossings
    for p3 in range(a2 + 1, a1):
        if p3 in sd.left_strands:
            q3 = sd.left_strands[p3]
            if b2 < q3 < b1:
                return Bimodule.Element()

    # next, check for black-orange double-crossings
    orange_strands_double_crossed = []

    # for each orange strand, we want to check if it double-crosses either of the black strands
    for orange in range(1, len(sd.etangle.signs)):
        left_y_pos = sd.etangle.left_y_pos(orange)
        if left_y_pos is None:
            continue
        middle_y_pos = sd.etangle.middle_y_pos(orange)
        times_crossed_p1 = 0
        times_crossed_p2 = 0

        # count how many times this orange strand crosses the p1 -> q2 black strand
        if (a1 < left_y_pos) ^ (min(a1, b1) < left_y_pos):
            times_crossed_p1 += 1
        if (min(a1, b1) < left_y_pos) ^ (b1 < left_y_pos):
            times_crossed_p1 += 1
        if (b1 < left_y_pos) ^ (b1 < middle_y_pos):
            times_crossed_p1 += 1

        # count how many times this orange strand crosses the p2 -> q1 black strand
        if (a2 < left_y_pos) ^ (min(a1, b1) < left_y_pos):
            times_crossed_p2 += 1
        if (min(a1, b1) < left_y_pos) ^ (b2 < left_y_pos):
            times_crossed_p2 += 1
        if (b2 < left_y_pos) ^ (b2 < middle_y_pos):
            times_crossed_p2 += 1

        # if either black strand is double-crossed, add this orange strand to the list
        if times_crossed_p1 > 1 or times_crossed_p2 > 1:
            orange_strands_double_crossed += [orange]

    # turn the list of orange strands into a coefficient
    for orange in orange_strands_double_crossed:
        if sd.etangle.middle_signs()[orange] == -1:
            c *= sd.etangle.strand_index_to_variable(orange)
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


def dmixed_case1(sd: StrandDiagram, b1: int, b2: int) -> Bimodule.Element:
    c = sd.etangle.polyring.one()
    a1 = sd.right_strands[b1]
    a2 = sd.right_strands[b2]
    sd_out = StrandDiagram(sd.etangle, sd.left_strands, swap_values(sd.right_strands, b1, b2))

    for b3 in sd.left_strands.values():
        if b2 < b3 < b1:
            return Bimodule.Element()

    for b3 in sd.right_strands.keys():
        a3 = sd.right_strands[b3]
        if (b3 < b1 and a3 > a1) or (b3 > b2 and a3 < a2):
            return Bimodule.Element()

    for orange in range(1, len(sd.etangle.signs)):
        if b2 < sd.etangle.middle_y_pos(orange) < b1:
            if sd.etangle.middle_signs()[orange] == 1:
                if sd.etangle.left_strand_straight(orange) or True:  # ???
                    return Bimodule.Element()
                elif sd.etangle.right_y_pos(orange) > a1:
                    c *= sd.etangle.strand_index_to_variable(orange)
                elif sd.etangle.right_y_pos(orange) < a2:
                    c *= sd.etangle.strand_index_to_variable(orange)
            else:
                if sd.etangle.left_strand_straight(orange) or True:  # ???
                    c *= sd.etangle.strand_index_to_variable(orange)
                elif sd.etangle.right_y_pos(orange) > a1:
                    return Bimodule.Element()
                elif sd.etangle.right_y_pos(orange) < a2:
                    return Bimodule.Element()

    return Bimodule.Element({sd_out: c})


def dmixed_case2(sd: StrandDiagram, b1: int, b2: int) -> Bimodule.Element:
    c = sd.etangle.polyring.one()
    a1 = sd.left_y_pos(b1)
    a2 = sd.left_y_pos(b2)
    sd_out = StrandDiagram(sd.etangle, swap_values(sd.left_strands, a1, a2), sd.right_strands)

    for b3 in sd.right_strands.keys():
        if b2 < b3 < b1:
            return Bimodule.Element()

    for b3 in sd.left_strands.values():
        a3 = sd.left_y_pos(b3)
        if b2 < b3 < b1 and (a3 > a2 or a3 < a1):
            return Bimodule.Element()

    for orange in range(1, len(sd.etangle.signs)):
        if b2 < sd.etangle.middle_y_pos(orange) < b1:
            if sd.etangle.middle_signs()[orange] == 1:
                if sd.etangle.right_strand_straight(orange) or True:  # ???
                    c *= sd.etangle.strand_index_to_variable(orange)
                elif sd.etangle.left_y_pos(orange) > a2:
                    return Bimodule.Element()
                elif sd.etangle.left_y_pos(orange) < a1:
                    return Bimodule.Element()
            else:
                if sd.etangle.right_strand_straight(orange) or True:  # ???
                    return Bimodule.Element()
                elif sd.etangle.left_y_pos(orange) > a2:
                    c *= sd.etangle.strand_index_to_variable(orange)
                elif sd.etangle.left_y_pos(orange) < a1:
                    c *= sd.etangle.strand_index_to_variable(orange)

    return Bimodule.Element({sd_out: c})


def dmixed_case3(sd: StrandDiagram, b1: int, b2: int) -> Bimodule.Element:
    c = sd.etangle.polyring.one()
    a1 = sd.left_y_pos(b1)
    a2 = sd.right_y_pos(b2)
    new_left_strands = dict(sd.left_strands)
    new_left_strands[a1] = b2
    new_right_strands = dict(sd.right_strands)
    del new_right_strands[b2]
    new_right_strands[b1] = a2
    sd_out = StrandDiagram(sd.etangle, new_left_strands, new_right_strands)

    for b3 in sd.left_strands.values():
        a3 = sd.left_y_pos(b3)
        if b2 < b3 < b1 and a3 < a1:
            return Bimodule.Element()

    for b3 in sd.right_strands.keys():
        a3 = sd.right_y_pos(b3)
        if b2 < b3 < b1 and a3 < a2:
            return Bimodule.Element()

    for orange in range(1, len(sd.etangle.signs)):
        if b2 < sd.etangle.middle_y_pos(orange) < b1:
            if sd.etangle.middle_signs()[orange] == 1:
                if sd.etangle.left_y_pos(orange) < a1:
                    return Bimodule.Element()
                elif sd.etangle.right_y_pos(orange) < a2:
                    c *= sd.etangle.strand_index_to_variable(orange)
            else:
                if sd.etangle.left_y_pos(orange) < a1:
                    c *= sd.etangle.strand_index_to_variable(orange)
                elif sd.etangle.right_y_pos(orange) < a2:
                    return Bimodule.Element()

    return Bimodule.Element({sd_out: c})


def dmixed_case4(sd: StrandDiagram, b1: int, b2: int) -> Bimodule.Element:
    c = sd.etangle.polyring.one()
    a1 = sd.right_y_pos(b1)
    a2 = sd.left_y_pos(b2)
    new_left_strands = dict(sd.left_strands)
    new_left_strands[a2] = b1
    new_right_strands = dict(sd.right_strands)
    del new_right_strands[b1]
    new_right_strands[b2] = a1
    sd_out = StrandDiagram(sd.etangle, new_left_strands, new_right_strands)

    for b3 in sd.left_strands.values():
        a3 = sd.left_y_pos(b3)
        if b2 < b3 < b1 and a3 > a2:
            return Bimodule.Element()

    for b3 in sd.right_strands.keys():
        a3 = sd.right_y_pos(b3)
        if b2 < b3 < b1 and a3 > a1:
            return Bimodule.Element()

    for orange in range(1, len(sd.etangle.signs)):
        if b2 < sd.etangle.middle_y_pos(orange) < b1:
            if sd.etangle.middle_signs()[orange] == 1:
                if sd.etangle.left_y_pos(orange) > a2:
                    return Bimodule.Element()
                elif sd.etangle.right_y_pos(orange) > a1:
                    c *= sd.etangle.strand_index_to_variable(orange)
            else:
                if sd.etangle.left_y_pos(orange) > a2:
                    c *= sd.etangle.strand_index_to_variable(orange)
                elif sd.etangle.right_y_pos(orange) > a1:
                    return Bimodule.Element()

    return Bimodule.Element({sd_out: c})


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
