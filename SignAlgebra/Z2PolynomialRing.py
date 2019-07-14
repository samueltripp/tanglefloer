from __future__ import annotations
from operator import add
from multimethod import multimethod
from typing import Set, Dict, FrozenSet, Tuple, overload, Union


class Z2PolynomialRing:
    def __init__(self, variables):
        self.variables = variables

    def zero(self) -> Z2Polynomial:
        return Z2Polynomial(self, set())

    def one(self) -> Z2Polynomial:
        return Z2Polynomial(self, {Z2Monomial(self, tuple([0 for j in range(len(self.variables))]))})

    def __getitem__(self, item) -> Z2Polynomial:
        i = self.variables.index(item)
        return Z2Monomial(self, tuple([1 if i == j else 0 for j in range(len(self.variables))])).to_polynomial()


# knows addition
# passes multiplication and equality off to Z2Monomial x Z2Monomial
class Z2Polynomial:
    @multimethod
    def __init__(self, ring: Z2PolynomialRing, terms: FrozenSet):
        self.ring = ring
        self.terms = terms

    @multimethod
    def __init__(self, ring: Z2PolynomialRing, terms: Set):
        self.ring = ring
        self.terms = frozenset(terms)

    def degree(self):
        if len(self.terms) == 0:
            return 0
        else:
            degrees = {term.degree() for term in self.terms}
            if len(degrees) > 1:
                raise Exception('non-homogeneous polynomial')
            else:
                return degrees.pop()

    def __add__(self, other: Z2Polynomial) -> Z2Polynomial:
        return Z2Polynomial(self.ring, self.terms ^ other.terms)

    @multimethod
    def __mul__(self, other):
        return other.__rmul__(self)

    @multimethod
    def __mul__(self, other: Z2Polynomial) -> Z2Polynomial:
        out = self.ring.zero()

        for term1 in self.terms:
            for term2 in other.terms:
                out += (term1 * term2).to_polynomial()

        return out

    @multimethod
    def __eq__(self, other: Z2Polynomial):
        return self.terms == other.terms

    def __hash__(self):
        return hash(self.terms)

    def __repr__(self):
        out = ''
        for term in self.terms:
            out += str(term) + ' + '
        return out[:-3] or str(0)


# knows multiplication, equality
# passes adding back to Z2Polynomial
class Z2Monomial:
    def __init__(self, ring, powers: Tuple):
        self.ring = ring
        self.powers = powers

    def degree(self):
        return sum(self.powers)

    def to_polynomial(self):
        return Z2Polynomial(self.ring, {self})

    def __mul__(self, other: Z2Monomial) -> Z2Monomial:
        return Z2Monomial(self.ring, tuple(map(add, self.powers, other.powers)))

    def __eq__(self, other: Z2Monomial):
        return self.ring == other.ring and self.powers == other.powers

    def __hash__(self):
        return hash(self.ring) + hash(self.powers)

    def __repr__(self):
        out = ''
        for (var, power) in enumerate(self.powers):
            if power == 1:
                out += self.ring.variables[var] + '.'
            elif power > 1:
                out += self.ring.variables[var] + '^' + str(power) + '.'
        return out[:-1] or str(1)
