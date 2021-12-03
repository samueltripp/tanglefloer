from __future__ import annotations
from typing import Tuple

from frozendict import frozendict
from multimethod import multimethod

from Functions.Functions import simplify_coefficients
from SignAlgebra.AMinus import AMinus
from SignAlgebra.Z2PolynomialRing import Z2Polynomial


class TensorAlgebra:
    def __init__(self, algebra: AMinus):
        self.algebra = algebra

    def zero(self):
        return TensorAlgebra.Element(self, {})

    def one_generator(self):
        return TensorAlgebra.Generator(self, tuple())

    def one(self):
        return self.one_generator().to_element()

    def __hash__(self):
        return hash(self.algebra)

    def __eq__(self, other: TensorAlgebra):
        return self.algebra == other.algebra

    class Element:
        def __init__(self, tensor_algebra, coefficients):
            self.tensor_algebra = tensor_algebra
            self.coefficients = frozendict(simplify_coefficients(coefficients))

        def to_algebra(self):
            return AMinus.Element(self.tensor_algebra.algebra,
                                  {g.to_algebra(): c for g, c in self.coefficients.items()})

        @multimethod
        def __add__(self, other: TensorAlgebra.Element) -> TensorAlgebra.Element:
            new_coefficients = dict(self.coefficients)
            for g in other.coefficients:
                if g in self.coefficients:
                    new_coefficients[g] = self.coefficients[g] + other.coefficients[g]
                else:
                    new_coefficients[g] = other.coefficients[g]
            return TensorAlgebra.Element(self.tensor_algebra, new_coefficients)

        @multimethod
        def __add__(self, other: TensorAlgebra.Generator) -> TensorAlgebra.Element:
            return self + other.to_element()

        # scalar multiplication
        def __rmul__(self, other: Z2Polynomial) -> TensorAlgebra.Element:
            assert other.ring == self.tensor_algebra.algebra.ring
            new_coefficients = {}
            for g, c in self.coefficients.items():
                new_coefficients[g] = other * c
            return TensorAlgebra.Element(self.tensor_algebra, new_coefficients)

        # the tensor product A^(x)i (x) A -> A^(x)i+1
        @multimethod
        def __pow__(self, other: AMinus.Element) -> TensorAlgebra.Element:
            out = self.tensor_algebra.zero()

            for g1, c1 in self.coefficients.items():
                for g2, c2 in other.coefficients.items():
                    if g1.right_idempotent() is None or g2.left_idempotent() is None \
                            or g1.right_idempotent() == g2.left_idempotent():
                        out += (c1 * c2) * (g1 ** g2)

            return out

        @multimethod
        def __pow__(self, other: AMinus.Generator) -> TensorAlgebra.Element:
            return self ** other.to_element()

        @multimethod
        def __pow__(self, other):
            return other.__rpow__(self)

        @multimethod
        def __rpow__(self, other: AMinus.Generator) -> TensorAlgebra.Element:
            return other.to_element() ** self

        # the tensor product A (x) A^(x)i -> A^(x)i+1
        @multimethod
        def __rpow__(self, other: AMinus.Element) -> TensorAlgebra.Element:
            out = self.tensor_algebra.zero()

            for g1, c1 in other.coefficients.items():
                for g2, c2 in self.coefficients.items():
                    if g1.right_idempotent() is None or g2.left_idempotent() is None \
                            or g1.right_idempotent() == g2.left_idempotent():
                        out += (c1 * c2) * (g1 ** g2)

            return out

        @multimethod
        def __eq__(self, other: TensorAlgebra.Element) -> bool:
            return self.tensor_algebra == other.tensor_algebra and self.coefficients == other.coefficients

        @multimethod
        def __eq__(self, other: TensorAlgebra.Generator) -> bool:
            return self == other.to_element()

        def __hash__(self):
            return hash((self.tensor_algebra, self.coefficients))

        def __repr__(self) -> str:
            return str(dict(self.coefficients))

    class Generator:
        def __init__(self, tensor_algebra, factors: Tuple):
            self.tensor_algebra = tensor_algebra
            self.factors = factors

        def to_element(self) -> TensorAlgebra.Element:
            return TensorAlgebra.Element(self.tensor_algebra, {self: self.tensor_algebra.algebra.ring.one()})

        def to_algebra(self):
            assert len(self.factors) == 1
            return self.factors[0]

        def num_factors(self):
            return len(self.factors)

        def left_idempotent(self) -> AMinus.Element | None:
            if len(self.factors) == 0:
                return None
            return self.factors[0].left_idempotent()

        def right_idempotent(self) -> AMinus.Element | None:
            if len(self.factors) == 0:
                return None
            return self.factors[-1].right_idempotent()

        def __add__(self, other) -> TensorAlgebra.Element:
            return self.to_element() + other

        def __radd__(self, other) -> TensorAlgebra.Element:
            return other + self.to_element()

        def __rmul__(self, other: Z2Polynomial) -> TensorAlgebra.Element:
            return other * self.to_element()

        # tensor product
        @multimethod
        def __pow__(self, other: AMinus.Generator) -> TensorAlgebra.Generator:
            assert self.right_idempotent() is None or other.left_idempotent() is None \
                or self.right_idempotent() == other.left_idempotent()
            return TensorAlgebra.Generator(self.tensor_algebra, self.factors + (other,))

        # tensor product
        @multimethod
        def __pow__(self, other: AMinus.Element) -> TensorAlgebra.Element:
            return self.to_element() ** other

        @multimethod
        def __pow__(self, other: TensorAlgebra.Generator) -> TensorAlgebra.Generator:
            assert self.tensor_algebra == other.tensor_algebra
            assert self.right_idempotent() is None or other.left_idempotent() is None \
                or self.right_idempotent() == other.left_idempotent()
            return TensorAlgebra.Generator(self.tensor_algebra, self.factors + other.factors)

        @multimethod
        def __pow__(self, other):
            return other.__rpow__(self)

        # tensor product
        @multimethod
        def __rpow__(self, other: AMinus.Generator) -> TensorAlgebra.Generator:
            assert self.right_idempotent() is None or other.left_idempotent() is None \
                    or self.right_idempotent() == other.left_idempotent()
            return TensorAlgebra.Generator(self.tensor_algebra, (other,) + self.factors)

        # tensor product
        @multimethod
        def __rpow__(self, other: AMinus.Element) -> TensorAlgebra.Element:
            return other ** self.to_element()

        def __str__(self):
            return str(self.factors)

        def __repr__(self):
            return str(self.factors)

        @multimethod
        def __eq__(self, other: TensorAlgebra.Generator):
            return self.tensor_algebra == other.tensor_algebra and self.factors == other.factors

        @multimethod
        def __eq__(self, other):
            return self.to_element() == other

        def __hash__(self):
            return hash((self.tensor_algebra, self.factors))
