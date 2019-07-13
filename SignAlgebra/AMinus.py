from __future__ import annotations
from multimethod import multimethod
import copy
from SignAlgebra.Z2PolynomialRing import *


class AMinus:
    def __init__(self, sign_sequence):
        # store the overarching sign sequence, the list of positive indices, and the polynomial ring
        self.ss = sign_sequence
        self.positives = list(i for i in range(len(sign_sequence)) if sign_sequence[i] == 1)
        self.polyring = self.initpolyring()

    # construct the polynomial ring based on the number of positives
    def initpolyring(self):
        np = len(self.positives)
        return Z2PolynomialRing(['U%s' % p for p in range(1, np + 1)])

    def zero(self):
        return AMinusElement(self, {})

    def gen(self, strands):
        return AMinusGen(self, strands).to_element()

    @multimethod
    def diff(self, gen: AMinusGen) -> AMinusElement:
        return self.gen_diff(gen)

    @multimethod
    def diff(self, elt: AMinusElement) -> AMinusElement:
        # sum the differentials of all the strand diagrams
        out = self.zero()
        for gen, coefficient in elt.coefficients.items():
            out += coefficient * self.gen_diff(gen)
        return out

    def gen_mult(self, gen1, gen2):
        strands1 = gen1.strands
        strands2 = gen2.strands
        # check if the ends don't match
        if set(strands1.values()) != strands2.keys():
            return self.zero()

        # check if black strands double cross
        for key1 in strands1.keys():
            for key2 in strands1.keys():
                if (strands2[strands1[key1]] > strands2[strands1[key2]] and strands1[key1] < strands1[key2]) \
                        or (strands2[strands1[key1]] < strands2[strands1[key2]] and strands1[key1] > strands1[key2]):
                    return self.zero()

        # count double crossed orange strands
        c = self.polyring.one()
        for key in strands1.keys():
            if strands1[key] > key:
                checkrange = range(max(key, strands2[strands1[key]]), strands1[key])
            else:
                checkrange = range(strands1[key], min(key, strands2[strands1[key]]))

            for i in checkrange:
                if i not in self.positives:
                    return self.zero()
                else:
                    c *= self.polyring['U' + str(self.positives.index(i) + 1)]

        # construct the new generator
        strands = {}
        for key in strands1.keys():
            strands[key] = strands2[strands1[key]]

        return c * AMinusGen(self, strands)

    def gen_diff(self, gen: AMinusGen) -> AMinusElement:
        # find all strands that cross, and resolve them
        out = self.zero()
        for s1, t1 in gen.strands.items():
            for s2, t2 in gen.strands.items():
                if s2 < s1 and t2 > t1:
                    out += AMinusElement(self, self.resolve(gen, s1, s2))
        return out

    def resolve(self, gen: AMinusGen, i, j):
        # check if black strands double cross
        for s in gen.strands.keys() & set(range(j, i)):
            if gen.strands[i] < gen.strands[s] < gen.strands[j]:
                return {}

        # construct the output generator
        out = AMinusGen(self, copy.deepcopy(gen.strands))
        out.strands[i] = gen.strands[j]
        out.strands[j] = gen.strands[i]

        # calculate the appropriate coefficient
        checkrange = range(max(gen.strands[i], j), min(i, gen.strands[j]))
        c = self.polyring.one()
        for i in checkrange:
            if i not in self.positives:
                return {}
            else:
                c = c * self.polyring['U' + str(self.positives.index(i) + 1)]

        return {out: c}

    @multimethod
    def twoalexander(self, elt: AMinusGen):
        return self.twoalex_gen(elt, self.polyring.one())

    @multimethod
    def twoalexander(self, elt: AMinusElement):
        firstkey = set(elt.coefficients.keys()).pop()
        twoalex = self.twoalex_gen(firstkey, elt.coefficients[firstkey])
        for key in elt.coefficients.keys():
            if self.twoalex_gen(key, elt.coefficients[key]) != twoalex:
                return None
        return twoalex

    def twoalex_gen(self, gen: AMinusGen, coeff):
        strands = dict(gen.strands)
        twoalex = -2 * coeff.degree()
        for key in strands.keys():
            if key < strands[key]:
                checkrange = range(key, strands[key])
            else:
                checkrange = range(strands[key], key)
            for k in checkrange:
                twoalex = twoalex + (-1) * self.ss[k]

        return twoalex

    def maslov(self, elt: AMinusElement):
        firstkey = set(elt.coefficients.keys()).pop
        maslov = self.maslov_gen(firstkey, elt.coefficients[firstkey])
        for key in elt.coefficients.keys():
            if self.maslov_gen(key, elt.coefficients[key]) != maslov:
                return None
        return maslov

    def maslov_gen(self, gen, coeff):
        strands = gen.strands
        maslov = -2 * coeff.degree()
        for key in strands.keys():
            for keyb in strands.keys():
                if keyb < key and strands[keyb] > strands[key]:
                    maslov = maslov + 1
            if key < strands[key]:
                checkrange = range(key, strands[key])
            else:
                checkrange = range(strands[key], key)
            for k in checkrange:
                if k in self.positives:
                    maslov = maslov - 1

        return maslov


class AMinusElement:
    def __init__(self, algebra, coefficients):
        self.algebra = algebra
        self.coefficients = {gen: coefficient for (gen, coefficient) in coefficients.items() if coefficient != 0}

    def __add__(self, other: AMinusElement):
        out_coefficients = {}
        for gen in self.coefficients.keys() - other.coefficients.keys():
            out_coefficients[gen] = self.coefficients[gen]
        for gen in self.coefficients.keys() & other.coefficients.keys():
            out_coefficients[gen] = self.coefficients[gen] + other.coefficients[gen]
        for gen in other.coefficients.keys() - self.coefficients.keys():
            out_coefficients[gen] = other.coefficients[gen]
        return AMinusElement(self.algebra, out_coefficients)

    @multimethod
    def __mul__(self, other: int) -> AMinusElement:
        out = self.algebra.zero()
        for (gen, coefficient) in self.coefficients.items():
            out += (other * coefficient) * gen
        return out

    @multimethod
    def __mul__(self, other: Z2Polynomial) -> AMinusElement:
        out = self.algebra.zero()
        for (gen, coefficient) in self.coefficients.items():
            out += (other * coefficient) * gen
        return out

    @multimethod
    def __mul__(self, other: AMinusElement):
        out = self.algebra.zero()
        for (gen1, coefficient1) in self.coefficients.items():
            for (gen2, coefficient2) in other.coefficients.items():
                out += (coefficient1 * coefficient2) * (gen1 * gen2)

        return out

    def __rmul__(self, other: Union[int, Z2Polynomial]) -> AMinusElement:
        return self * other

    def __eq__(self, other: AMinusElement):
        return self.algebra == other.algebra and self.coefficients == other.coefficients

    def __repr__(self):
        return str(self.coefficients)


# this only exists because I want to use strand diagrams as keys in algebra elements.
# strand diagrams are dictionaries, however, and you can't use dictionaries as keys in dictionaries.
# hence, wrap them in a class.
class AMinusGen:
    def __init__(self, algebra, strands):
        self.algebra = algebra
        self.strands = strands

    def to_element(self) -> AMinusElement:
        return AMinusElement(self.algebra, {self: 1})

    @multimethod
    def __mul__(self, other: int) -> AMinusElement:
        return AMinusElement(self.algebra, {self: other})

    @multimethod
    def __mul__(self, other: Z2Polynomial) -> AMinusElement:
        return AMinusElement(self.algebra, {self: other})

    @multimethod
    def __mul__(self, other: AMinusGen) -> AMinusElement:
        return self.algebra.gen_mult(self, other)

    @multimethod
    def __rmul__(self, other: int):
        return self * other

    @multimethod
    def __rmul__(self, other: Z2Polynomial):
        return self * other

    def __eq__(self, other):
        return self.algebra == other.algebra and self.strands == other.strands

    def __hash__(self):
        return hash(self.algebra) + hash(frozenset(self.strands.items()))

    def __repr__(self):
        return str(self.strands)
