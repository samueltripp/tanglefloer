from __future__ import annotations
from operator import add
from typing import Set, Dict, FrozenSet, Tuple, overload, Union


class Z2PolynomialRing:
    def __init__(self, variables):
        self.variables = variables

    def zero(self) -> Z2Polynomial:
        return Z2Polynomial(self, set())

    def one(self) -> Z2Monomial:
        return Z2Monomial(self, tuple([0 for j in range(len(self.variables))]))

    def __getitem__(self, item) -> Z2Monomial:
        i = self.variables.index(item)
        return Z2Monomial(self, tuple([1 if i == j else 0 for j in range(len(self.variables))]))


class Z2Polynomial:
    def __init__(self, ring, terms: Set[Z2Monomial]):
        self.ring = ring
        self.terms = terms

    def __add__(self, other: Union[Z2Polynomial, Z2Monomial]) -> Z2Polynomial:
        if isinstance(other, Z2Monomial):
            return self + Z2Polynomial(self.ring, {other})
        elif isinstance(other, Z2Polynomial):
            return Z2Polynomial(self.ring, self.terms.symmetric_difference(other.terms))

    def __mul__(self, other: Union[Z2Monomial, Z2Polynomial]) -> Z2Polynomial:
        if isinstance(other, Z2Monomial):
            return self * Z2Polynomial(self.ring, {other})
        elif isinstance(other, Z2Polynomial):
            out = self.ring.zero()

            for term1 in self.terms:
                for term2 in other.terms:
                    out += term1 * term2

            return out

    def __eq__(self, other):
        if isinstance(other, Z2Monomial):
            return self == Z2Polynomial(self.ring, {other})
        elif isinstance(other, Z2Polynomial):
            return self.terms == other.terms
        else:
            return False

    def __hash__(self):
        return hash(self.terms)

    def __repr__(self):
        out = ''
        for term in self.terms:
            out += str(term) + ' + '
        return out[:-3] or str(0)


class Z2Monomial:
    def __init__(self, ring, powers: Tuple):
        self.ring = ring
        self.powers = powers

    def degree(self):
        return sum(self.powers)

    def __add__(self, other: Union[Z2Monomial, Z2Polynomial]) -> Z2Polynomial:
        if isinstance(other, Z2Monomial):
            return Z2Polynomial(self.ring, {self}) + Z2Polynomial(self.ring, {other})
        elif isinstance(other, Z2Polynomial):
            return Z2Polynomial(self.ring, {self}) + other

    def __mul__(self, other):
        if isinstance(other, Z2Monomial):
            return Z2Monomial(self.ring, tuple(map(add, self.powers, other.powers)))
        elif isinstance(other, Z2Polynomial):
            return Z2Polynomial(self.ring, {self}) * other

    def __eq__(self, other):
        if isinstance(other, Z2Monomial):
            return self.ring == other.ring and self.powers == other.powers
        elif isinstance(other, Z2Polynomial):
            return Z2Polynomial(self.ring, {self}) == other
        else:
            return False

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
