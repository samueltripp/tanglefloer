from __future__ import annotations

from Tangles.Tangle import *
from Modules.Bimodule import *
from Modules.ETangleStrands import *


def type_da(etangle: ETangle) -> TypeDA:
    strands = [ETangleStrands(etangle, left_strands, right_strands)
               for left_strands, right_strands in
               enumerate_gens([etangle.left_points(), etangle.middle_points(), etangle.right_points()])]
    gens = [x.to_generator() for x in strands]
    maps = \
        sum((delta1_1(x) for x in strands), []) + \
        sum((delta1_2(x, a) for x in strands
             for a in etangle.right_algebra.left_gens(list(x.right_strands.values()))), [])

    return TypeDA(etangle.left_algebra, etangle.right_algebra, gens, maps)


def delta1_1(x: ETangleStrands) -> List[Bimodule.Edge]:
    out = []
    out += [Bimodule.Edge(x.to_generator(), y.to_generator(), c, (x.left_idempotent(),), tuple()) for y, c in d_plus(x).d.items()]
    out += [Bimodule.Edge(x.to_generator(), y.to_generator(), c, (x.left_idempotent(),), tuple()) for y, c in d_minus(x).d.items()]
    out += [Bimodule.Edge(x.to_generator(), y.to_generator(), c, (x.left_idempotent(),), tuple()) for y, c in d_mixed(x).d.items()]
    out += delta_ell(x)
    return out


def delta1_2(x: ETangleStrands, a: AMinus.Element) -> List[Bimodule.Edge]:
    b = m2(x, a)
    out = []
    for y, c in b.d.items():
        out += [Bimodule.Edge(x.to_generator(), y.to_generator(), c, (x.left_idempotent(),), (a,))]
    return out


def d_plus(x: ETangleStrands) -> Bimodule.Element:
    out = Bimodule.Element()
    for black1, black2 in itertools.combinations(x.right_strands.keys(), 2):
        b1 = min(black1, black2)
        b2 = max(black1, black2)
        if x.right_y_pos(b1) > x.right_y_pos(b2):
            # smooth the crossing and add it to the output
            out += smooth_right_crossing(x, b1, b2)
    return out


def d_minus(x: ETangleStrands) -> Bimodule.Element:
    out = Bimodule.Element()
    for black1, black2 in itertools.combinations(x.left_strands.values(), 2):
        b1 = min(black1, black2)
        b2 = max(black1, black2)
        # if the black strands don't cross
        if x.left_y_pos(b1) < x.left_y_pos(b2):
            # introduce a crossing and add it to the output
            out += introduce_left_crossing(x, b1, b2)
    return out


def d_mixed(x: ETangleStrands) -> Bimodule.Element:
    out = Bimodule.Element()

    for b1 in x.right_strands.keys():
        for b2 in x.right_strands.keys():
            # if two black strands don't cross on the right
            if b1 < b2 and x.right_y_pos(b1) < x.right_y_pos(b2):
                out += d_mixed_case_1(x, b1, b2)

    for b1 in x.left_strands.values():
        for b2 in x.left_strands.values():
            # if two black strands cross on the left
            if b1 < b2 and x.left_y_pos(b1) > x.left_y_pos(b2):
                out += d_mixed_case_2(x, b1, b2)

    for b1 in x.right_strands.keys():
        for b2 in x.left_strands.values():
            if b1 < b2:
                out += d_mixed_case_3(x, b1, b2)

    for b1 in x.left_strands.values():
        for b2 in x.right_strands.keys():
            if b1 < b2:
                out += d_mixed_case_4(x, b1, b2)

    return out


def delta_ell(x: ETangleStrands) -> List[Bimodule.Edge]:
    out = []

    unoccupied = set(x.etangle.left_points()) - set(x.left_strands.keys())

    for a1 in unoccupied:
        for a2 in unoccupied:
            if a1 < a2:
                e = delta_ell_case_1(x, a1, a2)
                if e:
                    out += [e]

    for a1 in x.left_strands.keys():
        for a2 in x.left_strands.keys():
            if a1 < a2 and x.left_strands[a1] > x.left_strands[a2]:
                e = delta_ell_case_2(x, a1, a2)
                if e:
                    out += [e]

    for a1 in unoccupied:
        for a2 in x.left_strands.keys():
            if a1 < a2:
                e = delta_ell_case_3(x, a1, a2)
                if e:
                    out += [e]

    for a1 in x.left_strands.keys():
        for a2 in unoccupied:
            if a1 < a2:
                e = delta_ell_case_4(x, a1, a2)
                if e:
                    out += [e]

    return out


@multimethod
def m2(x: ETangleStrands, a: AMinus.Element) -> Bimodule.Element:
    out = Bimodule.Element()
    for gen, coefficient in a.coefficients.items():
        out += x.etangle.from_right_algebra(gen.algebra, coefficient) * m2(x, gen)
    return out


@multimethod
def m2(x: ETangleStrands, a: AMinus.Gen) -> Bimodule.Element:
    # if the sign sequences do not match, return 0
    if x.etangle.right_signs() != a.algebra.ss:
        return Bimodule.Element()

    # if the strands cannot be merged, return 0
    if set(x.right_strands.values()) != set(a.strands.keys()):
        return Bimodule.Element()

    orange_strands = {}
    orange_signs = {}
    for orange in range(1, len(x.etangle.signs)):
        if x.etangle.right_y_pos(orange):
            orange_strands[orange] = (
                x.etangle.middle_y_pos(orange), x.etangle.right_y_pos(orange), x.etangle.right_y_pos(orange)
            )
        orange_signs[orange] = x.etangle.middle_signs()[orange]
    black_strands = {}
    for black in x.right_strands.keys():
        black_strands[black] = (black, x.right_strands[black], a.strands[x.right_strands[black]])

    c = x.etangle.polyring.one()
    sd = StrandDiagram(orange_strands, orange_signs, black_strands)
    powers = sd.figure_6_relations()
    if powers is None:
        return Bimodule.Element()
    for orange, power in powers.items():
        for _ in range(power):
            c *= x.etangle.strand_index_to_variable(orange)

    new_right_strands = {sd.black_left_pos(black): sd.black_right_pos(black) for black in x.right_strands.keys()}
    x_out = ETangleStrands(x.etangle, x.left_strands, new_right_strands)
    return Bimodule.Element({x_out: c})


def smooth_right_crossing(x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
    a1 = x.right_y_pos(b1)
    a2 = x.right_y_pos(b2)

    orange_strands = {}
    orange_signs = {}
    for orange in range(1, len(x.etangle.signs)):
        if x.etangle.right_y_pos(orange):
            orange_strands[orange] = (
                x.etangle.middle_y_pos(orange), x.etangle.middle_y_pos(orange), x.etangle.right_y_pos(orange)
            )
        orange_signs[orange] = x.etangle.middle_signs()[orange]
    black_strands = {}
    for black in x.right_strands.keys():
        if black == b1:
            black_strands[b1] = (b1, min(a1, b2) + .1, a2)
        elif black == b2:
            black_strands[b2] = (b2, min(a1, b2) + .2, a1)
        else:
            black_strands[black] = (black, black, x.right_y_pos(black))

    c = x.etangle.polyring.one()
    sd = StrandDiagram(orange_strands, orange_signs, black_strands)
    powers = sd.figure_6_relations()
    if powers is None:
        return Bimodule.Element()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power

    new_right_strands = swap_values(x.right_strands, b1, b2)
    x_out = ETangleStrands(x.etangle, x.left_strands, new_right_strands)
    return Bimodule.Element({x_out: c})


def introduce_left_crossing(x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
    a1 = x.left_y_pos(b1)
    a2 = x.left_y_pos(b2)

    orange_strands = {}
    orange_signs = {}
    for orange in range(1, len(x.etangle.signs)):
        if x.etangle.left_y_pos(orange):
            orange_strands[orange] = (
                x.etangle.left_y_pos(orange), x.etangle.middle_y_pos(orange), x.etangle.middle_y_pos(orange)
            )
        orange_signs[orange] = x.etangle.middle_signs()[orange]
    black_strands = {}
    for black in x.left_strands.values():
        if black == b1:
            black_strands[b1] = (a1, min(a2, b2) - .25, b1)
        elif black == b2:
            black_strands[b2] = (a2, min(a2, b2) + .25, b2)
        else:
            black_strands[black] = (x.left_y_pos(black), black, black)

    c = x.etangle.polyring.one()
    sd = StrandDiagram(orange_strands, orange_signs, black_strands)
    powers = sd.figure_7_relations()
    if powers is None:
        return Bimodule.Element()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power

    new_left_strands = swap_values(x.left_strands, a1, a2)
    x_out = ETangleStrands(x.etangle, new_left_strands, x.right_strands)
    return Bimodule.Element({x_out: c})


# swap the values associated to the given keys
def swap_values(d: Dict, k1, k2) -> Dict:
    d_out = dict(d)
    v1 = d[k1]
    v2 = d[k2]
    d_out[k1] = v2
    d_out[k2] = v1
    return d_out


def d_mixed_case_1(x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.to_strand_diagram().figure_8_case_1b(b1, b2)
    if powers is None:
        return Bimodule.Element()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    sd_out = ETangleStrands(x.etangle, x.left_strands, swap_values(x.right_strands, b1, b2))
    return Bimodule.Element({sd_out: c})


def d_mixed_case_2(x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.to_strand_diagram().figure_8_case_2b(b1, b2)
    if powers is None:
        return Bimodule.Element()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    a1 = x.left_y_pos(b1)
    a2 = x.left_y_pos(b2)
    sd_out = ETangleStrands(x.etangle, swap_values(x.left_strands, a1, a2), x.right_strands)
    return Bimodule.Element({sd_out: c})


def d_mixed_case_3(x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.to_strand_diagram().figure_8_case_3b(b1, b2)
    if powers is None:
        return Bimodule.Element()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    a1 = x.right_y_pos(b1)
    a2 = x.left_y_pos(b2)
    new_left_strands = dict(x.left_strands)
    new_left_strands[a2] = b1
    new_right_strands = dict(x.right_strands)
    del new_right_strands[b1]
    new_right_strands[b2] = a1
    sd_out = ETangleStrands(x.etangle, new_left_strands, new_right_strands)
    return Bimodule.Element({sd_out: c})


def d_mixed_case_4(x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.to_strand_diagram().figure_8_case_4b(b1, b2)
    if powers is None:
        return Bimodule.Element()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    a1 = x.left_y_pos(b1)
    a2 = x.right_y_pos(b2)
    new_left_strands = dict(x.left_strands)
    new_left_strands[a1] = b2
    new_right_strands = dict(x.right_strands)
    del new_right_strands[b2]
    new_right_strands[b1] = a2
    sd_out = ETangleStrands(x.etangle, new_left_strands, new_right_strands)
    return Bimodule.Element({sd_out: c})


def delta_ell_case_1(x: ETangleStrands, a1: int, a2: int) -> Optional[Bimodule.Edge]:
    c = x.etangle.polyring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_1a(a1, a2)
    if powers is None:
        return None
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    b1 = a1
    b2 = a2
    idempotent_strands = x.left_idempotent_strands()
    elt_strands = swap_values(idempotent_strands, b1, b2)
    elt = x.etangle.left_algebra.gen(elt_strands)
    return Bimodule.Edge(x.to_generator(), x.to_generator(), c, (elt,), tuple())


def delta_ell_case_2(x: ETangleStrands, a1: int, a2: int) -> Optional[Bimodule.Edge]:
    c = x.etangle.polyring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_2a(a1, a2)
    if powers is None:
        return None
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    elt = x.left_idempotent()
    y = ETangleStrands(x.etangle, swap_values(x.left_strands, a1, a2), x.right_strands)
    return Bimodule.Edge(x.to_generator(), y.to_generator(), c, (elt,), tuple())


def delta_ell_case_3(x: ETangleStrands, a1: int, a2: int) -> Optional[Bimodule.Edge]:
    c = x.etangle.polyring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_3a(a1, a2)
    if powers is None:
        return None
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    b1 = a1
    b2 = x.left_strands[a2]
    elt_strands = x.left_idempotent_strands()
    elt_strands[b1] = a2
    new_left_strands = dict(x.left_strands)
    del new_left_strands[a2]
    new_left_strands[a1] = b2
    elt = x.etangle.left_algebra.gen(elt_strands)
    y = ETangleStrands(x.etangle, new_left_strands, x.right_strands)
    return Bimodule.Edge(x.to_generator(), y.to_generator(), c, (elt,), tuple())


def delta_ell_case_4(x: ETangleStrands, a1: int, a2: int) -> Optional[Bimodule.Edge]:
    c = x.etangle.polyring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_4a(a1, a2)
    if powers is None:
        return None
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    b1 = x.left_strands[a1]
    b2 = a2
    elt_strands = x.left_idempotent_strands()
    elt_strands[b2] = a1
    new_left_strands = dict(x.left_strands)
    del new_left_strands[a1]
    new_left_strands[a2] = b1
    elt = x.etangle.left_algebra.gen(elt_strands)
    y = ETangleStrands(x.etangle, new_left_strands, x.right_strands)
    return Bimodule.Edge(x.to_generator(), y.to_generator(), c, (elt,), tuple())


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
