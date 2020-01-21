from __future__ import annotations

from abc import ABC, abstractmethod
from collections import defaultdict
from functools import lru_cache

from heapdict import heapdict
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
from Functions.Functions import simplify_coefficients
from pathos.pools import ProcessPool


# Base class for Type D, A, DD, AA, DA, and AD structures
class Module(ABC):
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
        self.reducible_edges = heapdict()

        for x in self.graph:
            for y in self.graph[x]:
                if self.edge_is_reducible(x, y):
                    self.reducible_edges[(x, y)] = self.reduction_cost(x, y)

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
        self.reducible_edges = heapdict()

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

        if self.edge_is_reducible(x, y):
            self.reducible_edges[(x, y)] = self.reduction_cost(x, y)
        elif (x, y) in self.reducible_edges:
            del self.reducible_edges[(x, y)]

    @abstractmethod
    def edge_is_reducible(self, x, y) -> bool:
        pass

    def reduction_cost(self, x, y) -> int:
        return len(self.graph.in_edges(y)) * len(self.graph.out_edges(x))

    def reduced(self) -> Module:
        components = self.decomposed()
        for component in components:
            component.reduce_component()
        return self.direct_sum(components)

    def pool_reduced(self) -> Module:
        components = self.decomposed()
        pool = ProcessPool(12)
        components_reduced = pool.map(lambda m: m.pool_reduce_component(), components)
        for component in components_reduced:
            component.restore_graph()
        return self.direct_sum(components_reduced)

    # O(V^2 + E)
    def restore_graph(self) -> None:
        self.graph.add_nodes_from(self.nodes)
        self.graph.add_edges_from(self.edges)

        for x in self.graph:
            for y in self.graph[x]:
                if self.edge_is_reducible(x, y):
                    self.reducible_edges[(x, y)] = self.reduction_cost(x, y)

    def pool_reduce_component(self):
        self.restore_graph()
        reducible_edge = self.get_reducible_edge()
        while reducible_edge is not None:
            self.reduce_edge(*reducible_edge)
            reducible_edge = self.get_reducible_edge()
        return self

    def reduce_component(self) -> Module:
        reducible_edge = self.get_reducible_edge()
        while reducible_edge is not None:
            self.reduce_edge(*reducible_edge)
            reducible_edge = self.get_reducible_edge()
        return self

    def get_reducible_edge(self):
        if len(self.reducible_edges) == 0:
            return None
        (x, y), _ = self.reducible_edges.popitem()
        while x not in self.graph or y not in self.graph:
            if len(self.reducible_edges) == 0:
                return None
            (x, y), _ = self.reducible_edges.popitem()
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
    def direct_sum(modules: List) -> Module:
        pass

    # add the given generator to this module
    def add_generator(self, generator: Module.TensorGenerator) -> None:
        self.graph.add_node(generator)

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
            if other.ring == self.module.ring:
                new_coefficients = {}
                for g, c in self.coefficients.items():
                    new_coefficients[g] = other * c
                return Module.TensorElement(self.module, self.i, self.j, new_coefficients)
            elif other.ring == self.module.left_algebra.ring:
                out = self.module.zero(self.i, self.j)
                for g, c in self.coefficients.items():
                    out += c * (other * g)
                return out

        # the tensor product A (x) (A^(x)i (x) M (x) A^(x)j) -> A^(x)i+1 (x) M (x) A^(x)j
        def __rpow__(self, other: AMinus.Element) -> Module.TensorElement:
            out = self.module.zero(self.i + 1, self.j)

            for g1, c1 in other.coefficients.items():
                for g2, c2 in self.coefficients.items():
                    if g1.right_idempotent() == g2.leftmost_idempotent():
                        if self.module.left_scalar_action is not None:
                            out += (self.module.left_scalar_action.apply(c1) * c2) * (g1 ** g2)
                        else:
                            out += c2 * (c1 * (g1 ** g2))

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
                            out += c1 * ((g1 ** g2) * c2)

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
                     left: Tuple[AMinus.Generator, ...] = None, left_monomial: Z2Monomial = None,
                     right: Tuple[AMinus.Generator, ...] = None, right_monomial: Z2Monomial = None):
            self.module = module
            self.key = key
            self.left_idempotent = left_idempotent
            self.right_idempotent = right_idempotent
            self.left = left or tuple()  # the tuple of generators of A^(x)i
            self.i = len(self.left)
            self.left_monomial = left_monomial or Z2Monomial(self.module.left_algebra.ring, {})
            self.right = right or tuple()  # the tuple of generators of A^(x)j
            self.j = len(self.right)
            self.right_monomial = right_monomial or Z2Monomial(self.module.right_algebra.ring, {})

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

        @multimethod
        def __rmul__(self, other: Z2Polynomial) -> Module.TensorElement:
            if other.ring == self.left_monomial.ring:
                out = self.module.zero(self.i, self.j)
                for term in other.terms:
                    out += (term * self).to_element()
                return out
            elif other.ring == self.module.ring:
                return other * self.to_element()

        @multimethod
        def __rmul__(self, other: Z2Monomial) -> Module.TensorGenerator:
            return Module.TensorGenerator(self.module, self.key, self.left_idempotent, self.right_idempotent,
                                          self.left, self.left_monomial * other,
                                          self.right, self.right_monomial)

        # tensor product
        @multimethod
        def __pow__(self, other: AMinus.Generator) -> Module.TensorGenerator:
            assert self.rightmost_idempotent() == other.left_idempotent()
            return Module.TensorGenerator(self.module, self.key, self.left_idempotent, self.right_idempotent,
                                          self.left, self.left_monomial,
                                          self.right + (other,), self.right_monomial)

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
            return str((self.left_monomial, self.left, self.key, self.right_monomial, self.right))

        def __repr__(self):
            return str((self.left_monomial, self.left, self.key, self.right_monomial, self.right))

        @multimethod
        def __eq__(self, other: Module.TensorGenerator):
            return self.module == other.module and \
                   self.key == other.key and \
                   self.left_idempotent == other.left_idempotent and \
                   self.right_idempotent == other.right_idempotent and \
                   self.left == other.left and \
                   self.left_monomial == other.left_monomial and \
                   self.right == other.right and \
                   self.right_monomial == other.right_monomial

        @multimethod
        def __eq__(self, other):
            return self.to_element() == other

        def __hash__(self):
            return hash((self.left_monomial, self.left, self.module, self.key,
                         self.left_idempotent, self.right_idempotent, self.right_monomial, self.right))
