from __future__ import annotations

from Modules.Bimodule import Bimodule
from Tangles.Tangle import *
from Modules.Bimodule import *
from Modules.ETangleStrands import *


def type_da(etangle: ETangle) -> TypeDA:
    out = TypeDA(etangle.polyring, etangle.left_algebra, etangle.right_algebra,
                 etangle.left_scalar_action, etangle.right_scalar_action)

    strands = [ETangleStrands(etangle, left_strands, right_strands)
               for left_strands, right_strands in
               enumerate_gens([etangle.left_points(), etangle.middle_points(), etangle.right_points()])]

    for x in strands:
        out.add_generator(x.to_generator(out))

    for x in strands:
        out.add_structure_map(x.to_generator(out), delta1_1(out, x))

        for a in etangle.right_algebra.left_gens(list(x.right_strands.values())):
            out.add_structure_map(x.to_generator(out) * a, delta1_2(out, x, a))

    return out


def delta1_1(module: TypeDA, x: ETangleStrands) -> Bimodule.Element:
    return x.left_idempotent() * (d_plus(module, x) + d_minus(module, x) + d_mixed(module, x)) + delta_ell(module, x)


def delta1_2(module: TypeDA, x: ETangleStrands, a: AMinus.Generator) -> Bimodule.Element:
    return x.left_idempotent() * m2(module, x, a)


def d_plus(module: Bimodule, x: ETangleStrands) -> Bimodule.Element:
    out = module.zero()
    for black1, black2 in itertools.combinations(x.right_strands.keys(), 2):
        b1 = min(black1, black2)
        b2 = max(black1, black2)
        if x.right_y_pos(b1) > x.right_y_pos(b2):
            # smooth the crossing and add it to the output
            out += smooth_right_crossing(module, x, b1, b2)
    return out


def d_minus(module: Bimodule, x: ETangleStrands) -> Bimodule.Element:
    out = module.zero()
    for black1, black2 in itertools.combinations(x.left_strands.values(), 2):
        b1 = min(black1, black2)
        b2 = max(black1, black2)
        # if the black strands don't cross
        if x.left_y_pos(b1) < x.left_y_pos(b2):
            # introduce a crossing and add it to the output
            out += introduce_left_crossing(module, x, b1, b2)
    return out


def d_mixed(module: Bimodule, x: ETangleStrands) -> Bimodule.Element:
    out = module.zero()

    for b1 in x.right_strands.keys():
        for b2 in x.right_strands.keys():
            # if two black strands don't cross on the right
            if b1 < b2 and x.right_y_pos(b1) < x.right_y_pos(b2):
                out += d_mixed_case_1(module, x, b1, b2)

    for b1 in x.left_strands.values():
        for b2 in x.left_strands.values():
            # if two black strands cross on the left
            if b1 < b2 and x.left_y_pos(b1) > x.left_y_pos(b2):
                out += d_mixed_case_2(module, x, b1, b2)

    for b1 in x.right_strands.keys():
        for b2 in x.left_strands.values():
            if b1 < b2:
                out += d_mixed_case_3(module, x, b1, b2)

    for b1 in x.left_strands.values():
        for b2 in x.right_strands.keys():
            if b1 < b2:
                out += d_mixed_case_4(module, x, b1, b2)

    return out


def delta_ell(module: Bimodule, x: ETangleStrands) -> Bimodule.Element:
    out = module.zero(1, 0)
    unoccupied = set(x.etangle.left_points()) - set(x.left_strands.keys())

    for a1 in unoccupied:
        for a2 in unoccupied:
            if a1 < a2:
                out += delta_ell_case_1(module, x, a1, a2)

    for a1 in x.left_strands.keys():
        for a2 in x.left_strands.keys():
            if a1 < a2 and x.left_strands[a1] > x.left_strands[a2]:
                out += delta_ell_case_2(module, x, a1, a2)

    for a1 in unoccupied:
        for a2 in x.left_strands.keys():
            if a1 < a2:
                out += delta_ell_case_3(module, x, a1, a2)

    for a1 in x.left_strands.keys():
        for a2 in unoccupied:
            if a1 < a2:
                out += delta_ell_case_4(module, x, a1, a2)

    return out


def m2(module: Bimodule, x: ETangleStrands, a: AMinus.Generator) -> Bimodule.Element:
    # if the sign sequences do not match, return 0
    if x.etangle.right_signs() != a.algebra.ss:
        return module.zero()

    # if the strands cannot be merged, return 0
    if set(x.right_strands.values()) != set(a.strands.keys()):
        return module.zero()

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
        return module.zero()
    for orange, power in powers.items():
        for _ in range(power):
            c *= x.etangle.strand_index_to_variable(orange)

    new_right_strands = {sd.black_left_pos(black): sd.black_right_pos(black) for black in x.right_strands.keys()}
    x_out = ETangleStrands(x.etangle, x.left_strands, new_right_strands)
    return c * x_out.to_generator(module)


def smooth_right_crossing(module: Bimodule, x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
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
        return module.zero()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power

    new_right_strands = swap_values(x.right_strands, b1, b2)
    x_out = ETangleStrands(x.etangle, x.left_strands, new_right_strands)
    return c * x_out.to_generator(module)


def introduce_left_crossing(module: Bimodule, x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
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
        return module.zero()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power

    new_left_strands = swap_values(x.left_strands, a1, a2)
    x_out = ETangleStrands(x.etangle, new_left_strands, x.right_strands)
    return c * x_out.to_generator(module)


# swap the values associated to the given keys
def swap_values(d: Dict, k1, k2) -> Dict:
    d_out = dict(d)
    v1 = d[k1]
    v2 = d[k2]
    d_out[k1] = v2
    d_out[k2] = v1
    return d_out


def d_mixed_case_1(module: Bimodule, x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.to_strand_diagram().figure_8_case_1b(b1, b2)
    if powers is None:
        return module.zero()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    x_out = ETangleStrands(x.etangle, x.left_strands, swap_values(x.right_strands, b1, b2))
    return c * x_out.to_generator(module)


def d_mixed_case_2(module: Bimodule, x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.to_strand_diagram().figure_8_case_2b(b1, b2)
    if powers is None:
        return module.zero()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    a1 = x.left_y_pos(b1)
    a2 = x.left_y_pos(b2)
    x_out = ETangleStrands(x.etangle, swap_values(x.left_strands, a1, a2), x.right_strands)
    return c * x_out.to_generator(module)


def d_mixed_case_3(module: Bimodule, x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.to_strand_diagram().figure_8_case_3b(b1, b2)
    if powers is None:
        return module.zero()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    a1 = x.right_y_pos(b1)
    a2 = x.left_y_pos(b2)
    new_left_strands = dict(x.left_strands)
    new_left_strands[a2] = b1
    new_right_strands = dict(x.right_strands)
    del new_right_strands[b1]
    new_right_strands[b2] = a1
    x_out = ETangleStrands(x.etangle, new_left_strands, new_right_strands)
    return c * x_out.to_generator(module)


def d_mixed_case_4(module: Bimodule, x: ETangleStrands, b1: int, b2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.to_strand_diagram().figure_8_case_4b(b1, b2)
    if powers is None:
        return module.zero()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    a1 = x.left_y_pos(b1)
    a2 = x.right_y_pos(b2)
    new_left_strands = dict(x.left_strands)
    new_left_strands[a1] = b2
    new_right_strands = dict(x.right_strands)
    del new_right_strands[b2]
    new_right_strands[b1] = a2
    x_out = ETangleStrands(x.etangle, new_left_strands, new_right_strands)
    return c * x_out.to_generator(module)


def delta_ell_case_1(module: Bimodule, x: ETangleStrands, a1: int, a2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_1a(a1, a2)
    if powers is None:
        return module.zero(1, 0)
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    b1 = a1
    b2 = a2
    idempotent_strands = x.left_idempotent_strands()
    elt_strands = swap_values(idempotent_strands, b1, b2)
    elt = x.etangle.left_algebra.generator(elt_strands)
    return elt * (c * x.to_generator(module))


def delta_ell_case_2(module: Bimodule, x: ETangleStrands, a1: int, a2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_2a(a1, a2)
    if powers is None:
        return module.zero(1, 0)
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    elt = x.left_idempotent()
    y = ETangleStrands(x.etangle, swap_values(x.left_strands, a1, a2), x.right_strands)
    return elt * (c * y.to_generator(module))


def delta_ell_case_3(module: Bimodule, x: ETangleStrands, a1: int, a2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_3a(a1, a2)
    if powers is None:
        return module.zero(1, 0)
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    b1 = a1
    b2 = x.left_strands[a2]
    elt_strands = x.left_idempotent_strands()
    elt_strands[b1] = a2
    new_left_strands = dict(x.left_strands)
    del new_left_strands[a2]
    new_left_strands[a1] = b2
    elt = x.etangle.left_algebra.generator(elt_strands)
    y = ETangleStrands(x.etangle, new_left_strands, x.right_strands)
    return elt * (c * y.to_generator(module))


def delta_ell_case_4(module: Bimodule, x: ETangleStrands, a1: int, a2: int) -> Bimodule.Element:
    c = x.etangle.polyring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_4a(a1, a2)
    if powers is None:
        return module.zero(1, 0)
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    b1 = x.left_strands[a1]
    b2 = a2
    elt_strands = x.left_idempotent_strands()
    elt_strands[b2] = a1
    new_left_strands = dict(x.left_strands)
    del new_left_strands[a1]
    new_left_strands[a2] = b1
    elt = x.etangle.left_algebra.generator(elt_strands)
    y = ETangleStrands(x.etangle, new_left_strands, x.right_strands)
    return elt * (c * y.to_generator(module))


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
