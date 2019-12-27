from __future__ import annotations

from collections import defaultdict
from functools import lru_cache

from networkx import MultiDiGraph
import networkx as nx
from typing import Iterable
from pygraphviz import AGraph
from Modules import ETangleStrands
from SignAlgebra.AMinus import AMinus
from Modules.CTMinus import *
from multimethod import *
from frozendict import *
from SignAlgebra.Z2PolynomialRing import *


# Base class for Type D, A, DD, AA, DA, and AD structures
class Module:
    def __init__(self, ring: Z2PolynomialRing, left_algebra: Optional[AMinus], right_algebra: Optional[AMinus],
                 left_scalar_action: Optional[Z2PolynomialRing.Map],
                 right_scalar_action: Optional[Z2PolynomialRing.Map],
                 graph: MultiDiGraph = None):
        self.ring = ring
        self.left_algebra = left_algebra
        self.right_algebra = right_algebra
        self.left_scalar_action = left_scalar_action
        self.right_scalar_action = right_scalar_action
        self.graph = graph or MultiDiGraph()

    def __repr__(self) -> str:
        return str(self.__dict__)

    # reduce this bimodule using the cancellation lemma
    def reduced(self) -> Module:
        out = self
        reducible_edge = out.reducible_edge()
        while reducible_edge is not None:
            out = out.reduce_edge(*reducible_edge)
            reducible_edge = out.reducible_edge()
        return out

    def reducible_edge(self):
        pass

    def reduce_edge(self, x, y, k, d):
        pass

    # add the given generator to this module
    def add_generator(self, generator: Module.Generator) -> None:
        self.graph.add_node(generator)

    # returns the zero element of A^(x)i (x) M (x) A^(x)j
    def zero(self, i=0, j=0) -> Module.Element:
        return Module.Element(self, i, j, {})

    # represents an element of A^(x)i (x) M (x) A^(x)j
    class Element:
        # coefficients - {TensorGenerator: Z2Polynomial}
        def __init__(self, module: Module, i, j, coefficients=None):
            self.module = module
            self.i = i
            self.j = j
            if coefficients is None:
                coefficients = {}
            self.coefficients = frozendict(coefficients)
            self.simplify()

        # remove any generators with a coefficient of 0 in this element
        # note: modifies this element instead of returning a new one
        def simplify(self) -> None:
            new_coefficients = dict(self.coefficients)
            for g, c in self.coefficients.items():
                if c == c.ring.zero():
                    del new_coefficients[g]
            self.coefficients = frozendict(new_coefficients)

        # addition in A^(x)i (x) M (x) A^(x)j
        def __add__(self, other: Module.Element) -> Module.Element:
            assert self.i == other.i and self.j == other.j
            new_coefficients = dict(self.coefficients)
            for g in other.coefficients:
                if g in self.coefficients:
                    new_coefficients[g] = self.coefficients[g] + other.coefficients[g]
                else:
                    new_coefficients[g] = other.coefficients[g]
            return Module.Element(self.module, self.i, self.j, new_coefficients)

        # scalar multiplication in A^(x)i (x) M (x) A^(x)j
        def __rmul__(self, other: Z2Polynomial) -> Module.Element:
            new_coefficients = dict(self.coefficients)
            for g in new_coefficients:
                new_coefficients[g] = other * new_coefficients[g]
            return Module.Element(self.module, self.i, self.j, new_coefficients)

        # the tensor product A (x) (A^(x)i (x) M (x) A^(x)j) -> A^(x)i+1 (x) M (x) A^(x)j
        def __rpow__(self, other: AMinus.Element) -> Module.Element:
            out = self.module.zero(self.i + 1, self.j)

            for g1, c1 in other.coefficients.items():
                for g2, c2 in self.coefficients.items():
                    if g1.right_idempotent() == g2.leftmost_idempotent():
                        out += (self.module.left_scalar_action.apply(c1) * c2) * (g1 ** g2)

            return out

        # the tensor product (A^(x)i (x) M (x) A^(x)j) (x) A -> A^(x)i (x) M (x) A^(x)j+1
        def __pow__(self, other: AMinus.Element) -> Module.Element:
            out = self.module.zero(self.i, self.j + 1)

            for g1, c1 in self.coefficients.items():
                for g2, c2 in other.coefficients.items():
                    if g1.rightmost_idempotent() == g2.left_idempotent():
                        out += (c1 * self.module.right_scalar_action.apply(c2)) * (g1 ** g2)

            return out

        @multimethod
        def __eq__(self, other: Module.Element) -> bool:
            return self.module == other.module and self.coefficients == other.coefficients

        @multimethod
        def __eq__(self, other: Module.Generator) -> bool:
            return self == other.to_element()

        def __hash__(self):
            return hash((self.module, self.coefficients))

        def __repr__(self) -> str:
            return str(self.coefficients)

    # represents a generator of A^(x)i (x) M (x) A^(x)j
    class Generator:
        def __init__(self, module: Module, key, left_idempotent: AMinus.Generator, right_idempotent: AMinus.Generator,
                     left: Tuple[AMinus.Generator, ...] = None, right: Tuple[AMinus.Generator, ...] = None):
            self.module = module  # the module this generator belongs to
            self.key = key  # a unique id for the module element m
            self.left_idempotent = left_idempotent  # the left idempotent of m
            self.right_idempotent = right_idempotent  # the right idempotent of m
            self.left = left or tuple()  # the tuple of generators of A^(x)i
            self.right = right or tuple()  # the tuple of generators of A^(x)j

        # converts this generator to an actual element
        def to_element(self) -> Module.Element:
            return Module.Element(self.module, len(self.left), len(self.right), {self: self.module.ring.one()})

        # returns m
        def get_module_generator(self):
            return Module.Generator(self.module, self.key,
                                    self.left_idempotent, self.right_idempotent, tuple(), tuple())

        def leftmost_idempotent(self) -> AMinus.Generator:
            if len(self.left) == 0:
                return self.left_idempotent
            else:
                return self.left[0].left_idempotent()

        def rightmost_idempotent(self) -> AMinus.Generator:
            if len(self.right) == 0:
                return self.right_idempotent
            else:
                return self.right[-1].right_idempotent()

        # addition
        def __add__(self, other) -> Module.Element:
            return self.to_element() + other

        # tensor product
        @multimethod
        def __pow__(self, other: AMinus.Generator) -> Module.Generator:
            assert self.rightmost_idempotent() == other.left_idempotent()
            return Module.Generator(self.module, self.key, self.left_idempotent, self.right_idempotent,
                                    self.left, self.right + (other,))

        # tensor product
        @multimethod
        def __pow__(self, other: AMinus.Element) -> Module.Element:
            return self.to_element() ** other

        # tensor product
        # assumes other is a tuple of generators or elements
        @multimethod
        def __pow__(self, other) -> Module.Generator:
            out = self
            for gen in other:
                out = out ** gen
            return out

        # tensor product
        @multimethod
        def __rpow__(self, other: AMinus.Generator) -> Module.Generator:
            assert other.right_idempotent() == self.leftmost_idempotent()
            return Module.Generator(self.module, self.key, self.left_idempotent, self.right_idempotent,
                                    (other,) + self.left, self.right)

        # tensor product
        @multimethod
        def __rpow__(self, other: AMinus.Element) -> Module.Element:
            return other ** self.to_element()

        # tensor product
        # assumes other is a tuple of generators or elements
        @multimethod
        def __rpow__(self, other) -> Module.Generator:
            out = self
            for gen in reversed(other):
                out = gen ** out
            return out

        # scalar multiplication
        def __rmul__(self, other: Z2Polynomial) -> Module.Element:
            return other * self.to_element()

        def __str__(self):
            return str((self.left, self.key, self.right))

        def __repr__(self):
            return str((self.left, self.key, self.right))

        @multimethod
        def __eq__(self, other: Module.Generator):
            return self.module == other.module and \
                   self.key == other.key and \
                   self.left == other.left and \
                   self.right == other.right

        @multimethod
        def __eq__(self, other):
            return self.to_element() == other

        def __hash__(self):
            return hash((self.module, self.key, self.left, self.right))
