from __future__ import annotations
from typing import List

from Modules.StrandDiagram import StrandDiagram
from SignAlgebra.Z2PolynomialRing import *
from Functions.Functions import *


# represents the algebra A^-(P) for some sign sequence P
class AMinus:
    def __init__(self, sign_sequence):
        # the sign sequence
        self.ss = tuple(sign_sequence) if sign_sequence[0] is None else (None,) + tuple(sign_sequence)
        # the list of positive indices, which are important for variable indices
        self.positives = (None,) + tuple([i for i, s in enumerate(self.ss) if s is not None and s > 0])
        # the polynomial ring acting on this algebra
        self.polyring = Z2PolynomialRing(['U%s' % p for p in range(1, len(self.positives) + 1)])

    # the zero element in A^-(P)
    def zero(self):
        return AMinus.Element(self, {})

    # a convenient way to construct elements of A^-(P)
    def generator(self, strands) -> AMinus.Generator:
        return AMinus.Generator(self, strands)

    # returns all generators with the given points occupied on the left
    def left_gens(self, points) -> List[AMinus.Generator]:
        return [AMinus.Generator(self, inj) for inj in injections(points, list(range(0, len(self.ss))))]

    # returns all generators with the given points occupied on the right
    def right_gens(self, points) -> List[AMinus.Generator]:
        return [AMinus.Generator(self, invert_injection(inj))
                for inj in injections(points, list(range(1, len(self.ss))))]

    # returns the idempotent with the given points occupied
    def idempotent(self, points) -> AMinus.Generator:
        return AMinus.Generator(self, {p: p for p in points})

    def __eq__(self, other: AMinus) -> bool:
        return self.ss == other.ss

    def __hash__(self):
        return hash(self.ss)

    # represents an element of A^-(P)
    class Element:
        def __init__(self, algebra, coefficients):
            self.algebra = algebra
            self.coefficients = frozendict({gen: coefficient for (gen, coefficient)
                                            in coefficients.items() if coefficient != algebra.polyring.zero()})

        # addition
        def __add__(self, other: AMinus.Element):
            out_coefficients = {}
            for gen in self.coefficients.keys() - other.coefficients.keys():
                out_coefficients[gen] = self.coefficients[gen]
            for gen in self.coefficients.keys() & other.coefficients.keys():
                out_coefficients[gen] = self.coefficients[gen] + other.coefficients[gen]
            for gen in other.coefficients.keys() - self.coefficients.keys():
                out_coefficients[gen] = other.coefficients[gen]
            return AMinus.Element(self.algebra, out_coefficients)

        # the algebra multiplication
        @multimethod
        def __mul__(self, other: AMinus.Element) -> AMinus.Element:
            out = self.algebra.zero()
            for (gen1, coefficient1) in self.coefficients.items():
                for (gen2, coefficient2) in other.coefficients.items():
                    out += (coefficient1 * coefficient2) * (gen1 * gen2)

            return out

        @multimethod
        def __mul__(self, other: AMinus.Generator) -> AMinus.Element:
            return self * other.to_element()

        # the scalar multiplication
        def __rmul__(self, other: Z2Polynomial) -> AMinus.Element:
            out = self.algebra.zero()
            for (gen, coefficient) in self.coefficients.items():
                out += (other * coefficient) * gen
            return out

        # the differential operation
        def diff(self) -> AMinus.Element:
            # sum the differentials of all the strand diagrams
            out = self.algebra.zero()
            for gen, coefficient in self.coefficients.items():
                out += coefficient * gen.diff()
            return out

        # twice the Alexander grading of this element
        def two_alexander(self):
            firstkey = set(self.coefficients.keys()).pop()
            out = firstkey.two_alexander(self.coefficients[firstkey])
            for key in self.coefficients.keys():
                if key.two_alexander(self.coefficients[key]) != out:
                    raise Exception('non-homogeneous element')
            return out

        # the Maslov grading of this element
        def maslov(self):
            firstkey = set(self.coefficients.keys()).pop()
            out = firstkey.maslov(self.coefficients[firstkey])
            for key in self.coefficients.keys():
                if key.maslov(self.coefficients[key]) != out:
                    raise Exception('non-homogeneous element')
            return out

        # returns True if this element is an idempotent
        def is_idempotent(self):
            return self == self * self

        @multimethod
        def __eq__(self, other: AMinus.Element):
            return self.algebra == other.algebra and self.coefficients == other.coefficients

        @multimethod
        def __eq__(self, other: AMinus.Generator):
            return self == other.to_element()

        def __hash__(self):
            return hash((self.algebra, frozendict(self.coefficients)))

        def __repr__(self):
            return str(self.coefficients)

    # this only exists because I want to use strand diagrams as keys in algebra elements.
    # strand diagrams are dictionaries, however, and you can't use dictionaries as keys in dictionaries.
    # hence, wrap them in a class.
    class Generator(Element):
        def __init__(self, algebra, strands):
            self.algebra = algebra
            self.strands = frozendict(strands)
            super().__init__(algebra, {self: algebra.polyring.one()})

        def to_element(self) -> AMinus.Element:
            return AMinus.Element(self.algebra, {self: self.algebra.polyring.one()})

        # the left idempotent of this element
        def left_idempotent(self) -> AMinus.Generator:
            return self.algebra.idempotent(self.strands.keys())

        # the right idempotent of this element
        def right_idempotent(self) -> AMinus.Generator:
            return self.algebra.idempotent(self.strands.values())

        # the algebra multiplication
        @multimethod
        def __mul__(self, other: AMinus.Generator) -> AMinus.Element:
            strands1 = self.strands
            strands2 = other.strands
            # check if the ends don't match
            if set(strands1.values()) != strands2.keys():
                return self.algebra.zero()

            orange_strands = {orange: 3 * (orange - 1 / 2,) for orange in range(1, len(self.algebra.ss))}
            orange_signs = {orange: self.algebra.ss[orange] for orange in range(1, len(self.algebra.ss))}
            black_strands = {self.strands[black]: (black, self.strands[black], other.strands[self.strands[black]])
                             for black in self.strands.keys()}

            sd = StrandDiagram(orange_strands, orange_signs, black_strands)
            powers = sd.figure_6_relations()
            if powers is None:
                return self.algebra.zero()
            c = self.algebra.polyring.one()
            for orange, power in powers.items():
                if orange in self.algebra.positives:
                    p = self.algebra.positives.index(orange)
                    c *= self.algebra.polyring['U' + str(p)] ** power

            # construct the new generator
            strands = {}
            for key in strands1.keys():
                strands[key] = strands2[strands1[key]]

            return c * AMinus.Generator(self.algebra, strands)

        @multimethod
        def __mul__(self, other: AMinus.Element):
            return self.to_element() * other

        # we don't know how to multiply generators by anything else
        @multimethod
        def __mul__(self, other):
            return other.__rmul__(self)

        # scalar multiplication
        @multimethod
        def __rmul__(self, other: Z2Polynomial) -> AMinus.Element:
            return AMinus.Element(self.algebra, {self: other})

        @multimethod
        def __rmul__(self, other: AMinus.Element) -> AMinus.Element:
            return other * self.to_element()

        @multimethod
        def __rmul__(self, other):
            raise NotImplementedError()

        # the differential
        def diff(self) -> AMinus.Element:
            # find all strands that cross, and resolve them
            out = self.algebra.zero()
            for s1, t1 in self.strands.items():
                for s2, t2 in self.strands.items():
                    if s1 < s2 and t1 > t2:
                        out += self.smooth_crossing(s1, s2)
            return out

        # computes a single summand of the differential
        def smooth_crossing(self, i, j) -> AMinus.Element:
            orange_strands = {orange: 3 * (orange - 1 / 2,) for orange in range(1, len(self.algebra.ss))}
            orange_signs = {orange: self.algebra.ss[orange] for orange in range(1, len(self.algebra.ss))}
            black_strands = {}

            for black in self.strands.keys():
                if black == i:
                    black_strands[i] = (i, min(j, self.strands[i]) - .25, self.strands[j])
                elif black == j:
                    black_strands[j] = (j, min(j, self.strands[i]) + .25, self.strands[i])
                else:
                    black_strands[black] = (black, black, self.strands[black])

            sd = StrandDiagram(orange_strands, orange_signs, black_strands)
            powers = sd.figure_6_relations()
            c = self.algebra.polyring.one()
            if powers is None:
                return self.algebra.zero()
            for orange, power in powers.items():
                if orange in self.algebra.positives:
                    p = self.algebra.positives.index(orange)
                    c *= self.algebra.polyring['U' + str(p)] ** power

            # construct the new generator
            new_strands = dict(self.strands)
            new_strands[i] = self.strands[j]
            new_strands[j] = self.strands[i]

            return c * AMinus.Generator(self.algebra, new_strands)

        # twice the ALexander grading of this generator
        def two_alexander(self, coeff):
            strands = dict(self.strands)
            out = -2 * coeff.degree()
            for key in strands.keys():
                if key < strands[key]:
                    checkrange = range(key, strands[key])
                else:
                    checkrange = range(strands[key], key)
                for k in checkrange:
                    out = out + (-1) * self.algebra.ss[k]

            return out

        # the Maslov grading of this generator
        def maslov(self, coeff):
            strands = self.strands
            out = -2 * coeff.degree()
            for key in strands.keys():
                for keyb in strands.keys():
                    if keyb < key and strands[keyb] > strands[key]:
                        out += 1
                if key < strands[key]:
                    checkrange = range(key, strands[key])
                else:
                    checkrange = range(strands[key], key)
                for k in checkrange:
                    if k in self.algebra.positives:
                        out -= 1

            return out

        @multimethod
        def __eq__(self, other: AMinus.Generator):
            return self.algebra == other.algebra and self.strands == other.strands

        @multimethod
        def __eq__(self, other):
            return AMinus.Element.__eq__(self, other)

        def __hash__(self):
            return hash((self.algebra, self.strands))

        def __repr__(self):
            return dict_to_sorted_string(self.strands)
