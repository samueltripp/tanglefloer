from __future__ import annotations

from abc import ABC, abstractmethod
from collections import defaultdict
from functools import lru_cache

from heapdict import heapdict
from networkx import MultiDiGraph
import networkx as nx
from typing import Iterable, List, Optional
from pygraphviz import AGraph
from Modules import ETangleStrands
from SignAlgebra.AMinus import AMinus
from Modules.CTMinus import *
from multimethod import *
from frozendict import *
from SignAlgebra.Z2PolynomialRing import *
from Functions.Functions import simplify_coefficients
from pathos.pools import ProcessPool


# Base class for Type D, A, DD, AA, DA, and AD structures
class Module(ABC):
    def __init__(self, ring: Z2PolynomialRing, left_algebra: Optional[AMinus], right_algebra: Optional[AMinus],
                 left_scalar_action: Optional[Z2PolynomialRing.Map],
                 right_scalar_action: Optional[Z2PolynomialRing.Map],
                 graph: MultiDiGraph = None,
                 gradings: dict = None):
        self.ring = ring
        self.left_algebra = left_algebra
        self.right_algebra = right_algebra
        self.left_scalar_action = left_scalar_action
        self.right_scalar_action = right_scalar_action
        self.graph = graph or MultiDiGraph()
        self.gradings = gradings or {}
        if left_scalar_action is not None:
            assert left_algebra.ring == left_scalar_action.source, "left scalar action has wrong source"
            assert ring == left_scalar_action.target, "left scalar action has wrong target"
        if right_scalar_action is not None:
            assert right_algebra.ring == right_scalar_action.source, "right scalar action has wrong source"
            assert ring == right_scalar_action.target, "right scalar action has wrong target"

        if graph is not None:
            if gradings is None:
                raise AssertionError("Incomplete grading information")
            assert set(graph.nodes).issubset(set(gradings.keys())), "Incomplete grading information"

    def __getstate__(self):
        return self.ring, self.left_algebra, self.right_algebra, self.left_scalar_action, self.right_scalar_action, \
               list(self.graph.nodes), list(self.graph.edges(keys=True, data=True))

    def __setstate__(self, state):
        self.ring = state[0]
        self.left_algebra = state[1]
        self.right_algebra = state[2]
        self.left_scalar_action = state[3]
        self.right_scalar_action = state[4]
        self.nodes = state[5]
        self.edges = state[6]
        self.graph = MultiDiGraph()
        self.edge_priorities = False
        self.reducible_edges = None

    def __repr__(self) -> str:
        return str(self.__dict__)

    def add_edge(self, x, y, k, c):
        current = self.graph.get_edge_data(x, y, key=k)
        if current is None:
            self.graph.add_edge(x, y, key=k, c=self.ring.zero())
            current = self.graph.get_edge_data(x, y, key=k)
        current['c'] += c
        if current['c'] == self.ring.zero():
            self.graph.remove_edge(x, y, key=k)

    @abstractmethod
    def edge_is_reducible(self, x, y) -> bool:
        pass

    def reduce(self) -> Module:
        components = self.decomposed()
        if len(components) == 0:
            return self
        for component in components:
            component.reduce_component()
        return self.direct_sum(components)

    def reduce_component(self) -> Module:
        reducible_edge = self.get_reducible_edge()
        while reducible_edge is not None:
            self.reduce_edge(*reducible_edge)
            reducible_edge = self.get_reducible_edge()
        return self

    def get_reducible_edge(self):
        for x in self.graph:
            for y in self.graph[x]:
                if self.edge_is_reducible(x, y):
                    k, d = list(self.graph[x][y].items())[0]
                    return x, y, k, d

    @abstractmethod
    def reduce_edge(self, x, y, k, d):
        pass

    @abstractmethod
    def decomposed(self) -> List[Module]:
        pass

    @staticmethod
    @abstractmethod
    def direct_sum(modules: List[Module]) -> Module:
        pass

    # add the given generator to this module
    def add_generator(self, generator: Module.TensorGenerator,grading:List[int]) -> None:
        self.graph.add_node(generator)
        self.gradings[generator]=grading
        self.gradings[generator] = grading

    # returns the zero element of A^(x)i (x) M (x) A^(x)j
    def zero(self, i=0, j=0) -> Module.TensorElement:
        return Module.TensorElement(self, i, j)

    # represents an element of A^(x)i (x) M (x) A^(x)j
    class TensorElement:
        # coefficients - {TensorGenerator: Z2Polynomial}
        def __init__(self, module: Module, i, j, coefficients=None):
            self.module = module
            self.i = i
            self.j = j
            if coefficients is None:
                coefficients = {}
            self.coefficients = frozendict(simplify_coefficients(coefficients))

        # addition in A^(x)i (x) M (x) A^(x)j
        def __add__(self, other: Module.TensorElement) -> Module.TensorElement:
            assert self.i == other.i and self.j == other.j
            new_coefficients = dict(self.coefficients)
            for g in other.coefficients:
                if g in self.coefficients:
                    new_coefficients[g] = self.coefficients[g] + other.coefficients[g]
                else:
                    new_coefficients[g] = other.coefficients[g]
            return Module.TensorElement(self.module, self.i, self.j, new_coefficients)

        # scalar multiplication in A^(x)i (x) M (x) A^(x)j as a module over the base ring of M
        def __rmul__(self, other: Z2Polynomial) -> Module.TensorElement:
            if (self.module.left_algebra is not None) and (other.ring == self.module.left_algebra.ring) \
                    and (self.module.left_scalar_action is not None):
                other = self.module.left_scalar_action.apply(other)
            elif (self.module.right_algebra is not None) and (other.ring == self.module.right_algebra.ring) \
                    and (self.module.right_scalar_action is not None):
                other = self.module.right_scalar_action.apply(other)
            if other.ring == self.module.ring:
                new_coefficients = {}
                for g, c in self.coefficients.items():
                    new_coefficients[g] = other * c
                return Module.TensorElement(self.module, self.i, self.j, new_coefficients)
            else:
                raise Exception('no action by the given polynomial')

        @multimethod
        def __rpow__(self, other: AMinus.Generator) -> Module.TensorElement:
            return other.to_element() ** self

        # the tensor product A (x) (A^(x)i (x) M (x) A^(x)j) -> A^(x)i+1 (x) M (x) A^(x)j
        @multimethod
        def __rpow__(self, other: AMinus.Element) -> Module.TensorElement:
            out = self.module.zero(self.i + 1, self.j)

            for g1, c1 in other.coefficients.items():
                for g2, c2 in self.coefficients.items():
                    if g1.right_idempotent() == g2.leftmost_idempotent():
                        if self.module.left_scalar_action is not None:
                            out += (self.module.left_scalar_action.apply(c1) * c2) * (g1 ** g2)
                        else:
                            raise Exception('no left scalar action')

            return out

        # the tensor product (A^(x)i (x) M (x) A^(x)j) (x) A -> A^(x)i (x) M (x) A^(x)j+1
        def __pow__(self, other: AMinus.Element) -> Module.TensorElement:
            out = self.module.zero(self.i, self.j + 1)

            for g1, c1 in self.coefficients.items():
                for g2, c2 in other.coefficients.items():
                    if g1.rightmost_idempotent() == g2.left_idempotent():
                        if self.module.right_scalar_action is not None:
                            out += (c1 * self.module.right_scalar_action.apply(c2)) * (g1 ** g2)
                        else:
                            raise Exception('no right scalar action')

            return out

        @multimethod
        def __eq__(self, other: Module.TensorElement) -> bool:
            return self.module == other.module and self.coefficients == other.coefficients

        @multimethod
        def __eq__(self, other: Module.TensorGenerator) -> bool:
            return self == other.to_element()

        def __hash__(self):
            return hash((self.module, self.coefficients))

        def __repr__(self) -> str:
            return str(dict(self.coefficients))

    # represents a generator of A^(x)i (x) M (x) A^(x)j as a module over the base ring of M
    class TensorGenerator:
        def __init__(self, module: Module, key,
                     left_idempotent: AMinus.Generator, right_idempotent: AMinus.Generator,
                     left: Tuple[AMinus.Generator, ...] = None,
                     right: Tuple[AMinus.Generator, ...] = None):
            self.module = module
            self.key = key
            self.left_idempotent = left_idempotent
            self.right_idempotent = right_idempotent
            self.left = left or tuple()  # the tuple of generators of A^(x)i
            self.i = len(self.left)
            self.right = right or tuple()  # the tuple of generators of A^(x)j
            self.j = len(self.right)

        # converts this generator to an actual element
        def to_element(self) -> Module.TensorElement:
            return Module.TensorElement(self.module, len(self.left), len(self.right), {self: self.module.ring.one()})

        def get_module_generator(self) -> Module.TensorGenerator:
            return Module.TensorGenerator(self.module, self.key, self.left_idempotent, self.right_idempotent)

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
        def __add__(self, other) -> Module.TensorElement:
            return self.to_element() + other

        def __rmul__(self, other: Z2Polynomial) -> Module.TensorElement:
            return other * self.to_element()

        # tensor product
        @multimethod
        def __pow__(self, other: AMinus.Generator) -> Module.TensorGenerator:
            assert self.rightmost_idempotent() == other.left_idempotent()
            return Module.TensorGenerator(self.module, self.key, self.left_idempotent, self.right_idempotent,
                                          self.left, self.right + (other,))

        # tensor product
        @multimethod
        def __pow__(self, other: AMinus.Element) -> Module.TensorElement:
            return self.to_element() ** other

        # tensor product
        # assumes other is a tuple of generators or elements
        @multimethod
        def __pow__(self, other) -> Module.TensorGenerator:
            out = self
            for gen in other:
                out = out ** gen
            return out

        # tensor product
        @multimethod
        def __rpow__(self, other: AMinus.Generator) -> Module.TensorGenerator:
            assert other.right_idempotent() == self.leftmost_idempotent()
            return Module.TensorGenerator(self.module, self.key,
                                          self.left_idempotent, self.right_idempotent, (other,) + self.left, self.right)

        # tensor product
        @multimethod
        def __rpow__(self, other: AMinus.Element) -> Module.TensorElement:
            return other ** self.to_element()

        # tensor product
        # assumes other is a tuple of generators or elements
        @multimethod
        def __rpow__(self, other) -> Module.TensorGenerator:
            out = self
            for gen in reversed(other):
                out = gen ** out
            return out

        def __str__(self):
            return str((self.left, self.key, self.right))

        def __repr__(self):
            return str((self.left, self.key, self.right))

        @multimethod
        def __eq__(self, other: Module.TensorGenerator):
            return self.module == other.module and \
                   self.key == other.key and \
                   self.left_idempotent == other.left_idempotent and \
                   self.right_idempotent == other.right_idempotent and \
                   self.left == other.left and \
                   self.right == other.right

        @multimethod
        def __eq__(self, other):
            return self.to_element() == other

        def __hash__(self):
            return hash((self.left, self.module, self.key,
                         self.left_idempotent, self.right_idempotent, self.right))
