from __future__ import annotations

from Modules.StrandDiagram import StrandDiagram
from SignAlgebra.Z2PolynomialRing import *
from Tangles.Functions import *


class AMinus:
    def __init__(self, sign_sequence):
        # store the overarching sign sequence, the list of positive indices, and the polynomial ring
        self.ss = tuple(sign_sequence) if sign_sequence[0] is None else (None,) + tuple(sign_sequence)
        self.positives = (None,) + tuple([i for i, s in enumerate(self.ss) if s is not None and s > 0])
        self.polyring = Z2PolynomialRing(['U%s' % p for p in range(1, len(self.positives) + 1)])

    def zero(self):
        return AMinus.Element(self, {})

    def gen(self, strands):
        return AMinus.Gen(self, strands).to_element()

    # returns all generators with the given points occupied on the left
    def left_gens(self, points):
        return [self.gen(inj) for inj in injections(points, list(range(0, len(self.ss))))]

    # returns all generators with the given points occupied on the right
    def right_gens(self, points):
        return [self.gen(invert_injection(inj)) for inj in injections(points, list(range(1, len(self.ss))))]

    # returns the idempotent with the given points occupied
    def idempotent(self, points):
        return self.gen({p: p for p in points})

    class Element:
        def __init__(self, algebra, coefficients):
            self.algebra = algebra
            self.coefficients = {gen: coefficient for (gen, coefficient)
                                 in coefficients.items() if coefficient != algebra.polyring.zero()}
    
        def __add__(self, other: AMinus.Element):
            out_coefficients = {}
            for gen in self.coefficients.keys() - other.coefficients.keys():
                out_coefficients[gen] = self.coefficients[gen]
            for gen in self.coefficients.keys() & other.coefficients.keys():
                out_coefficients[gen] = self.coefficients[gen] + other.coefficients[gen]
            for gen in other.coefficients.keys() - self.coefficients.keys():
                out_coefficients[gen] = other.coefficients[gen]
            return AMinus.Element(self.algebra, out_coefficients)
    
        @multimethod
        def __mul__(self, other: Z2Polynomial) -> AMinus.Element:
            out = self.algebra.zero()
            for (gen, coefficient) in self.coefficients.items():
                out += (other * coefficient) * gen
            return out
    
        @multimethod
        def __mul__(self, other: AMinus.Element) -> AMinus.Element:
            out = self.algebra.zero()
            for (gen1, coefficient1) in self.coefficients.items():
                for (gen2, coefficient2) in other.coefficients.items():
                    out += (coefficient1 * coefficient2) * (gen1 * gen2)
    
            return out
    
        def __rmul__(self, other: Z2Polynomial) -> AMinus.Element:
            return self * other
    
        def diff(self) -> AMinus.Element:
            # sum the differentials of all the strand diagrams
            out = self.algebra.zero()
            for gen, coefficient in self.coefficients.items():
                out += coefficient * gen.diff()
            return out
    
        def two_alexander(self):
            firstkey = set(self.coefficients.keys()).pop()
            out = firstkey.two_alexander(self.coefficients[firstkey])
            for key in self.coefficients.keys():
                if key.two_alexander(self.coefficients[key]) != out:
                    raise Exception('non-homogeneous element')
            return out
    
        def maslov(self):
            firstkey = set(self.coefficients.keys()).pop()
            out = firstkey.maslov(self.coefficients[firstkey])
            for key in self.coefficients.keys():
                if key.maslov(self.coefficients[key]) != out:
                    raise Exception('non-homogeneous element')
            return out
    
        def __eq__(self, other: AMinus.Element):
            return self.algebra == other.algebra and self.coefficients == other.coefficients
    
        def __repr__(self):
            return str(self.coefficients)

    # this only exists because I want to use strand diagrams as keys in algebra elements.
    # strand diagrams are dictionaries, however, and you can't use dictionaries as keys in dictionaries.
    # hence, wrap them in a class.
    class Gen:
        def __init__(self, algebra, strands):
            self.algebra = algebra
            self.strands = strands

        def to_element(self) -> AMinus.Element:
            return AMinus.Element(self.algebra, {self: self.algebra.polyring.one()})

        @multimethod
        def __mul__(self, other: Z2Polynomial) -> AMinus.Element:
            return AMinus.Element(self.algebra, {self: other})

        @multimethod
        def __mul__(self, other: AMinus.Gen) -> AMinus.Element:
            strands1 = self.strands
            strands2 = other.strands
            # check if the ends don't match
            if set(strands1.values()) != strands2.keys():
                return self.algebra.zero()

            orange_strands = {orange: 3*(orange - 1/2,) for orange in range(1, len(self.algebra.ss))}
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
                    c *= self.algebra.polyring['U'+str(p)] ** power

            # construct the new generator
            strands = {}
            for key in strands1.keys():
                strands[key] = strands2[strands1[key]]

            return c * AMinus.Gen(self.algebra, strands)

        def __rmul__(self, other: Z2Polynomial):
            return self * other

        def diff(self) -> AMinus.Element:
            # find all strands that cross, and resolve them
            out = self.algebra.zero()
            for s1, t1 in self.strands.items():
                for s2, t2 in self.strands.items():
                    if s1 < s2 and t1 > t2:
                        out += self.smooth_crossing(s1, s2)
            return out

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
                    c *= self.algebra.polyring['U'+str(p)] ** power

            # construct the new generator
            new_strands = dict(self.strands)
            new_strands[i] = self.strands[j]
            new_strands[j] = self.strands[i]

            return c * self.algebra.gen(new_strands)

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

        def __eq__(self, other):
            return self.algebra == other.algebra and self.strands == other.strands

        def __hash__(self):
            return hash(self.algebra) + hash(frozenset(self.strands.items()))

        def __repr__(self):
            return str(self.strands)
