from __future__ import annotations
from operator import add

from frozendict import frozendict
from multimethod import multimethod
from typing import Set, FrozenSet, Tuple, Dict, Iterable


class Z2PolynomialRing:
    def __init__(self, variables: Iterable):
        self.variables = set(variables)

    def zero(self) -> Z2Polynomial:
        return Z2Polynomial(self, set())

    def one(self) -> Z2Polynomial:
        return Z2Polynomial(self, {Z2Monomial(self, {})})

    def __getitem__(self, item) -> Z2Polynomial:
        assert item in self.variables
        return Z2Monomial(self, {item: 1}).to_polynomial()

    # returns the two maps
    # self -> self (x) other
    # other -> self (x) other
    def tensor_inclusions(self, other: Z2PolynomialRing):
        product = Z2PolynomialRing([v + 'a' for v in self.variables] + [v + 'b' for v in other.variables])

        in1 = Z2PolynomialRing.Map(self, product, {v: v+'a' for v in self.variables})
        in2 = Z2PolynomialRing.Map(other, product, {v: v+'b' for v in other.variables})

        return in1, in2

    def __eq__(self, other: Z2PolynomialRing):
        return self.variables == other.variables

    def __hash__(self):
        return hash(frozenset(self.variables))

    # represents a map between polynomial rings that sends some variables to other variables
    # very limited in scope
    class Map:
        # mapping: {source_variable_index: target_variable_index}
        def __init__(self, source: Z2PolynomialRing, target: Z2PolynomialRing, mapping: Dict):
            self.source = source
            self.target = target
            self.mapping = mapping

        def apply(self, x: Z2Polynomial):
            assert x.ring == self.source

            y = self.target.zero()

            for x_term in x.terms:
                y += Z2Monomial(self.target, {self.mapping[var]: power for var, power in x_term.powers.items()}).to_polynomial()

            return y


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

    def degree(self) -> int:
        if len(self.terms) == 0:
            return 0
        else:
            degrees = {term.degree() for term in self.terms}
            if len(degrees) > 1:
                raise Exception('non-homogeneous polynomial')
            else:
                return degrees.pop()

    def __add__(self, other: Z2Polynomial) -> Z2Polynomial:
        assert self.ring == other.ring

        return Z2Polynomial(self.ring, self.terms ^ other.terms)

    @multimethod
    def __mul__(self, other):
        return other.__rmul__(self)

    @multimethod
    def __mul__(self, other: Z2Polynomial) -> Z2Polynomial:
        assert self.ring == other.ring

        out = self.ring.zero()

        for term1 in self.terms:
            for term2 in other.terms:
                out += (term1 * term2).to_polynomial()

        return out

    def __pow__(self, power: int) -> Z2Polynomial:
        out = self.ring.one()
        for _ in range(power):
            out *= self
        return out

    @multimethod
    def __eq__(self, other: Z2Polynomial) -> bool:
        return self.terms == other.terms

    def __hash__(self):
        return hash(self.terms)

    def __repr__(self) -> str:
        out = ''
        for term in self.terms:
            out += str(term) + ' + '
        return out[:-3] or str(0)


# knows multiplication, equality
# passes adding back to Z2Polynomial
class Z2Monomial:
    def __init__(self, ring, powers: Dict):
        self.ring = ring
        self.powers = frozendict(powers)

    def degree(self) -> int:
        return sum(self.powers.values())

    def to_polynomial(self) -> Z2Polynomial:
        return Z2Polynomial(self.ring, {self})

    def __mul__(self, other: Z2Monomial) -> Z2Monomial:
        return Z2Monomial(self.ring, {var: self.powers.get(var, 0) + other.powers.get(var, 0)
                                      for var in self.powers.keys() | other.powers.keys()})

    def __eq__(self, other: Z2Monomial) -> bool:
        return self.ring == other.ring and self.powers == other.powers

    def __hash__(self):
        return hash(self.ring) + hash(self.powers)

    def __repr__(self) -> str:
        out = ''
        for var, power in self.powers.items():
            if power == 1:
                out += str(var) + '.'
            elif power > 1:
                out += str(var) + '^' + str(power) + '.'
        return out[:-1] or str(1)
