from __future__ import annotations

from Modules.TypeDA import TypeDA
from Modules.Module import Module
from Tangles.Tangle import *
from Modules.Module import *
from Modules.ETangleStrands import *


def type_da(etangle: ETangle) -> TypeDA:
    out = TypeDA(etangle.ring, etangle.left_algebra, etangle.right_algebra, etangle.right_scalar_action)

    strands = [ETangleStrands(etangle, left_strands, right_strands)
               for left_strands, right_strands in
               enumerate_gens([etangle.left_points(), etangle.middle_points(), etangle.right_points()])]

    for x in strands:
        out.add_generator(x.to_generator(out))

    for x in strands:
        out.add_structure_map(x.to_generator(out), delta1_1(out, x))

        for a in etangle.right_algebra.left_gens(list(x.right_strands.values())):
            out.add_structure_map(x.to_generator(out) ** a, delta1_2(out, x, a))

    return out


def delta1_1(module: TypeDA, x: ETangleStrands) -> Module.TensorElement:
    return x.left_idempotent().to_element() ** (d_plus(module, x) + d_minus(module, x) + d_mixed(module, x)) + delta_ell(module, x)


def delta1_2(module: TypeDA, x: ETangleStrands, a: AMinus.Generator) -> Module.TensorElement:
    return x.left_idempotent().to_element() ** m2(module, x, a)


def d_plus(module: Module, x: ETangleStrands) -> Module.TensorElement:
    out = module.zero()
    for black1, black2 in itertools.combinations(x.right_strands.keys(), 2):
        b1 = min(black1, black2)
        b2 = max(black1, black2)
        if x.right_y_pos(b1) > x.right_y_pos(b2):
            # smooth the crossing and add it to the output
            out += smooth_right_crossing(module, x, b1, b2)
    return out


def d_minus(module: Module, x: ETangleStrands) -> Module.TensorElement:
    out = module.zero()
    for black1, black2 in itertools.combinations(x.left_strands.values(), 2):
        b1 = min(black1, black2)
        b2 = max(black1, black2)
        # if the black strands don't cross
        if x.left_y_pos(b1) < x.left_y_pos(b2):
            # introduce a crossing and add it to the output
            out += introduce_left_crossing(module, x, b1, b2)
    return out


def d_mixed(module: Module, x: ETangleStrands) -> Module.TensorElement:
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


def delta_ell(module: Module, x: ETangleStrands) -> Module.TensorElement:
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


def m2(module: Module, x: ETangleStrands, a: AMinus.Generator) -> Module.TensorElement:
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

    c = x.etangle.ring.one()
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


@multimethod
def smooth_right_crossing(module: Module, x: ETangleStrands, b1: int, b2: int) -> Module.TensorElement:
    a1 = x.right_y_pos(b1)
    a2 = x.right_y_pos(b2)

    crossing_range_1 = range(a1, b1+1) if a1 < b1 else range(b1, a1+1)
    crossing_range_2 = range(a2, b2+1) if a2 < b2 else range(b2, a2+1)
    possible_crossings = list(set(crossing_range_1) & set(crossing_range_2))
    sd = smooth_right_crossing(x, b1, b2, possible_crossings[0])
    for crossing in possible_crossings[1:]:
        new_sd = smooth_right_crossing(x, b1, b2, crossing)
        if new_sd.num_orange_black_crossings() < sd.num_orange_black_crossings():
            sd = new_sd

    c = x.etangle.ring.one()
    powers = sd.figure_6_relations()
    if powers is None:
        return module.zero()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power

    new_right_strands = swap_values(x.right_strands, b1, b2)
    x_out = ETangleStrands(x.etangle, x.left_strands, new_right_strands)
    return c * x_out.to_generator(module)


@multimethod
def smooth_right_crossing(x: ETangleStrands, b1: int, b2: int, crossing: int) -> StrandDiagram:
    a1 = x.right_y_pos(b1)
    a2 = x.right_y_pos(b2)

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
        if black == b1:
            black_strands[b1] = (b1, crossing + .1, a2)
        elif black == b2:
            black_strands[b2] = (b2, crossing + .2, a1)
        else:
            black_strands[black] = (black, black, x.right_y_pos(black))

    return StrandDiagram(orange_strands, orange_signs, black_strands)


@multimethod
def introduce_left_crossing(module: Module, x: ETangleStrands, b1: int, b2: int) -> Module.TensorElement:
    a1 = x.left_y_pos(b1)
    a2 = x.left_y_pos(b2)

    crossing_range_1 = range(a1, b2+1) if a1 < b2 else range(b2, a1+1)
    crossing_range_2 = range(a2, b1+1) if a2 < b1 else range(b1, a2+1)
    possible_crossings = list(set(crossing_range_1) & set(crossing_range_2))
    sd = introduce_left_crossing(x, b1, b2, possible_crossings[0])
    for crossing in possible_crossings[1:]:
        new_sd = introduce_left_crossing(x, b1, b2, crossing)
        if new_sd.num_orange_black_crossings() < sd.num_orange_black_crossings():
            sd = new_sd

    c = x.etangle.ring.one()
    powers = sd.figure_7_relations()
    if powers is None:
        return module.zero()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power

    new_left_strands = swap_values(x.left_strands, a1, a2)
    x_out = ETangleStrands(x.etangle, new_left_strands, x.right_strands)
    return c * x_out.to_generator(module)


@multimethod
def introduce_left_crossing(x: ETangleStrands, b1: int, b2: int, crossing: int) -> StrandDiagram:
    a1 = x.left_y_pos(b1)
    a2 = x.left_y_pos(b2)

    orange_strands = {}
    orange_signs = {}
    for orange in range(1, len(x.etangle.signs)):
        if x.etangle.left_y_pos(orange):
            # (left, middle, middle) breaks {0: 0, 1: 2} on a cup with (-1, 1)
            orange_strands[orange] = (
                x.etangle.left_y_pos(orange), x.etangle.left_y_pos(orange), x.etangle.middle_y_pos(orange)
            )
        orange_signs[orange] = x.etangle.middle_signs()[orange]
    black_strands = {}
    for black in x.left_strands.values():
        if black == b1:
            black_strands[b1] = (a1, crossing + .1, b1)
        elif black == b2:
            black_strands[b2] = (a2, crossing + .2, b2)
        else:
            black_strands[black] = (x.left_y_pos(black), black, black)

    return StrandDiagram(orange_strands, orange_signs, black_strands)


def d_mixed_case_1(module: Module, x: ETangleStrands, b1: int, b2: int) -> Module.TensorElement:
    c = x.etangle.ring.one()
    powers = x.to_strand_diagram().figure_8_case_1b(b1, b2)
    if powers is None:
        return module.zero()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    x_out = ETangleStrands(x.etangle, x.left_strands, swap_values(x.right_strands, b1, b2))
    return c * x_out.to_generator(module)


def d_mixed_case_2(module: Module, x: ETangleStrands, b1: int, b2: int) -> Module.TensorElement:
    c = x.etangle.ring.one()
    powers = x.to_strand_diagram().figure_8_case_2b(b1, b2)
    if powers is None:
        return module.zero()
    for orange, power in powers.items():
        c *= x.etangle.strand_index_to_variable(orange) ** power
    a1 = x.left_y_pos(b1)
    a2 = x.left_y_pos(b2)
    x_out = ETangleStrands(x.etangle, swap_values(x.left_strands, a1, a2), x.right_strands)
    return c * x_out.to_generator(module)


def d_mixed_case_3(module: Module, x: ETangleStrands, b1: int, b2: int) -> Module.TensorElement:
    c = x.etangle.ring.one()
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


def d_mixed_case_4(module: Module, x: ETangleStrands, b1: int, b2: int) -> Module.TensorElement:
    c = x.etangle.ring.one()
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


def delta_ell_case_1(module: Module, x: ETangleStrands, a1: int, a2: int) -> Module.TensorElement:
    c1 = module.left_algebra.ring.one()
    c2 = module.ring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_1a(a1, a2)
    if powers is None:
        return module.zero(1, 0)
    for orange, power in powers.items():
        if x.etangle.left_signs()[orange] == 1:
            c1 *= x.etangle.left_algebra_strand_index_to_variable(orange) ** power
        else:
            c2 *= x.etangle.strand_index_to_variable(orange) ** power
    b1 = a1
    b2 = a2
    idempotent_strands = x.left_idempotent_strands()
    elt_strands = swap_values(idempotent_strands, b1, b2)
    elt = x.etangle.left_algebra.generator(elt_strands)
    return (c1 * elt) ** (c2 * x.to_generator(module))


def delta_ell_case_2(module: Module, x: ETangleStrands, a1: int, a2: int) -> Module.TensorElement:
    c1 = module.left_algebra.ring.one()
    c2 = x.etangle.ring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_2a(a1, a2)
    if powers is None:
        return module.zero(1, 0)
    for orange, power in powers.items():
        if x.etangle.left_signs()[orange] == 1:
            c1 *= x.etangle.left_algebra_strand_index_to_variable(orange) ** power
        else:
            c2 *= x.etangle.strand_index_to_variable(orange) ** power
    elt = x.left_idempotent()
    y = ETangleStrands(x.etangle, swap_values(x.left_strands, a1, a2), x.right_strands)
    return (c1 * elt) ** (c2 * y.to_generator(module))


def delta_ell_case_3(module: Module, x: ETangleStrands, a1: int, a2: int) -> Module.TensorElement:
    c1 = module.left_algebra.ring.one()
    c2 = x.etangle.ring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_3a(a1, a2)
    if powers is None:
        return module.zero(1, 0)
    for orange, power in powers.items():
        if x.etangle.left_signs()[orange] == 1:
            c1 *= x.etangle.left_algebra_strand_index_to_variable(orange) ** power
        else:
            c2 *= x.etangle.strand_index_to_variable(orange) ** power
    b1 = a1
    b2 = x.left_strands[a2]
    elt_strands = x.left_idempotent_strands()
    elt_strands[b1] = a2
    new_left_strands = dict(x.left_strands)
    del new_left_strands[a2]
    new_left_strands[a1] = b2
    elt = x.etangle.left_algebra.generator(elt_strands)
    y = ETangleStrands(x.etangle, new_left_strands, x.right_strands)
    return (c1 * elt) ** (c2 * y.to_generator(module))


def delta_ell_case_4(module: Module, x: ETangleStrands, a1: int, a2: int) -> Module.TensorElement:
    c1 = module.left_algebra.ring.one()
    c2 = x.etangle.ring.one()
    powers = x.idempotent_and_left_strands().figure_8_case_4a(a1, a2)
    if powers is None:
        return module.zero(1, 0)
    for orange, power in powers.items():
        if x.etangle.left_signs()[orange] == 1:
            c1 *= x.etangle.left_algebra_strand_index_to_variable(orange) ** power
        else:
            c2 *= x.etangle.strand_index_to_variable(orange) ** power
    b1 = x.left_strands[a1]
    b2 = a2
    elt_strands = x.left_idempotent_strands()
    elt_strands[b2] = a1
    new_left_strands = dict(x.left_strands)
    del new_left_strands[a1]
    new_left_strands[a2] = b1
    elt = x.etangle.left_algebra.generator(elt_strands)
    y = ETangleStrands(x.etangle, new_left_strands, x.right_strands)
    return (c1 * elt) ** (c2 * y.to_generator(module))


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
