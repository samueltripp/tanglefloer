from __future__ import annotations

from frozendict import frozendict
from multimethod import multimethod
from typing import Set, FrozenSet, Tuple, Dict, Iterable
from Functions.Functions import is_injection, invert_injection


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

    # represents a map between polynomial rings that sends some variables to other variables
    # very limited in scope
    class Map:
        # mapping: {source_variable_index: target_variable_index}
        def __init__(self, source: Z2PolynomialRing, target: Z2PolynomialRing, mapping: Dict):
            for s, t in mapping.items():
                assert s in source.variables and t in target.variables
            self.source = source
            self.target = target
            self.mapping = mapping
            self.retract = invert_injection(self.mapping) if is_injection(self.mapping) else None

        @staticmethod
        def identity(source: Z2PolynomialRing, target: Z2PolynomialRing) -> Z2PolynomialRing.Map:
            assert source.variables == target.variables

            return Z2PolynomialRing.Map(source, target, {v: v for v in source.variables})

        @multimethod
        def apply(self, x: Z2Polynomial) -> Z2Polynomial:
            assert x.ring == self.source

            y = self.target.zero()

            for x_term in x.terms:
                y += self.apply(x_term).to_polynomial()

            return y

        @multimethod
        def apply(self, x: Z2Monomial):
            assert x.ring == self.source

            powers = {}
            for var, power in x.powers.items():
                output_var = self.mapping[var]
                if output_var in powers:
                    powers[output_var] += power
                else:
                    powers[output_var] = power

            return Z2Monomial(self.target, powers)

        # applies f^{-1} to y if possible
        def retract(self, y: Z2Polynomial) -> Z2Polynomial:
            assert self.retract is not None and y.ring == self.target

            x = self.target.zero()

            for y_term in y.terms:
                x += Z2Monomial(self.source, {self.retract[var]: power
                                              for var, power in y_term.powers.items()}).to_polynomial()

            return x

        # returns the map x -> self.apply(other.apply(x))
        def compose(self, other: Z2PolynomialRing.Map) -> Z2PolynomialRing.Map:
            assert self.source == other.target
            return Z2PolynomialRing.Map(other.source, self.target,
                                        {var: self.mapping[other.mapping[var]] for var in other.mapping.keys()})

        # returns the two maps induced by the pushout of self and other
        def pushout_inclusions(self, other: Z2PolynomialRing.Map) -> Tuple:
            assert self.source == other.source
            left_vars = self.target.variables - set(self.mapping.values())
            common_vars = self.source.variables
            right_vars = other.target.variables - set(other.mapping.values())
            pushout = Z2PolynomialRing([v + 'a' for v in left_vars] +
                                       [v + 'b' for v in common_vars] +
                                       [v + 'c' for v in right_vars])

            in_left = Z2PolynomialRing.Map(self.target, pushout,
                                           {**{v: v + 'a' for v in left_vars},
                                            **{self.mapping[v]: v + 'b' for v in common_vars}})
            in_right = Z2PolynomialRing.Map(other.target, pushout,
                                            {**{v: v + 'c' for v in right_vars},
                                             **{other.mapping[v]: v + 'b' for v in common_vars}})

            return in_left, in_right


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

    @multimethod
    def __mul__(self, other):
        return other.__rmul__(self)

    @multimethod
    def __mul__(self, other: Z2Monomial) -> Z2Monomial:
        return Z2Monomial(self.ring, {var: self.powers.get(var, 0) + other.powers.get(var, 0)
                                      for var in self.powers.keys() | other.powers.keys()})

    @multimethod
    def __eq__(self, other: Z2Monomial) -> bool:
        return self.ring == other.ring and self.powers == other.powers

    @multimethod
    def __eq__(self, other) -> bool:
        return self.to_polynomial() == other

    def __hash__(self):
        return hash(self.ring) + hash(self.powers)

    def __repr__(self) -> str:
        out = ''
        for var, power in self.powers.items():
            if power == 1:
                out += str(var) + '*'
            elif power > 1:
                out += str(var) + '^' + str(power) + '.'
        return out[:-1] or str(1)
