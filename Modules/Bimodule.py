from __future__ import annotations

from collections import defaultdict
from functools import lru_cache

from networkx import MultiDiGraph
import networkx as nx
from typing import Iterable, Any
from pygraphviz import AGraph
from Modules import ETangleStrands
from SignAlgebra.AMinus import AMinus
from Modules.CTMinus import *
from multimethod import *
from frozendict import *
from SignAlgebra.Z2PolynomialRing import *


# Base class for Type DD, AA, DA, and AD structures
class Bimodule:
    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 left_scalar_action: Z2PolynomialRing.Map, right_scalar_action: Z2PolynomialRing.Map,
                 graph: MultiDiGraph=None):
        self.ring = ring
        self.left_algebra = left_algebra
        self.right_algebra = right_algebra
        self.left_scalar_action = left_scalar_action
        self.right_scalar_action = right_scalar_action
        self.graph = graph or MultiDiGraph()

    def __repr__(self) -> str:
        return str(self.__dict__)

    # returns the zero element of A^(x)i (x) M (x) A^(x)j
    def zero(self, i=0, j=0) -> Bimodule.Element:
        return Bimodule.Element(self, i, j, {})

    # represents an element of A^(x)i (x) M (x) A^(x)j
    class Element:
        # coefficients - {TensorGenerator: Z2Polynomial}
        def __init__(self, module: Bimodule, i, j, coefficients=None):
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
        def __add__(self, other: Bimodule.Element) -> Bimodule.Element:
            assert self.i == other.i and self.j == other.j
            new_coefficients = dict(self.coefficients)
            for g in other.coefficients:
                if g in self.coefficients:
                    new_coefficients[g] = self.coefficients[g] + other.coefficients[g]
                else:
                    new_coefficients[g] = other.coefficients[g]
            return Bimodule.Element(self.module, self.i, self.j, new_coefficients)

        # scalar multiplication in A^(x)i (x) M (x) A^(x)j
        def __rmul__(self, other: Z2Polynomial) -> Bimodule.Element:
            new_coefficients = dict(self.coefficients)
            for g in new_coefficients:
                new_coefficients[g] = other * new_coefficients[g]
            return Bimodule.Element(self.module, self.i, self.j, new_coefficients)

        # the tensor product A (x) (A^(x)i (x) M (x) A^(x)j) -> A^(x)i+1 (x) M (x) A^(x)j
        def __rpow__(self, other: AMinus.Element) -> Bimodule.Element:
            out = self.module.zero(self.i + 1, self.j)

            for g1, c1 in other.coefficients.items():
                for g2, c2 in self.coefficients.items():
                    if g1.right_idempotent() == g2.leftmost_idempotent():
                        out += (self.module.left_scalar_action.apply(c1) * c2) * (g1 ** g2)

            return out

        # the tensor product (A^(x)i (x) M (x) A^(x)j) (x) A -> A^(x)i (x) M (x) A^(x)j+1
        def __pow__(self, other: AMinus.Element) -> Bimodule.Element:
            out = self.module.zero(self.i, self.j + 1)

            for g1, c1 in self.coefficients.items():
                for g2, c2 in other.coefficients.items():
                    if g1.rightmost_idempotent() == g2.left_idempotent():
                        out += (c1 * self.module.right_scalar_action.apply(c2)) * (g1 ** g2)

            return out

        @multimethod
        def __eq__(self, other: Bimodule.Element) -> bool:
            return self.module == other.module and self.coefficients == other.coefficients

        @multimethod
        def __eq__(self, other: Bimodule.Generator) -> bool:
            return self == other.to_element()

        def __hash__(self):
            return hash((self.module, self.coefficients))

        def __repr__(self) -> str:
            return str(self.coefficients)

    # represents a generator of A^(x)i (x) M (x) A^(x)j
    class Generator:
        def __init__(self, module: Bimodule, key, left_idempotent: AMinus.Generator, right_idempotent: AMinus.Generator,
                     left: Tuple[AMinus.Generator, ...] = None, right: Tuple[AMinus.Generator, ...] = None):
            self.module = module  # the module this generator belongs to
            self.key = key  # a unique id for the module element m
            self.left_idempotent = left_idempotent  # the left idempotent of m
            self.right_idempotent = right_idempotent  # the right idempotent of m
            self.left = left or tuple()  # the tuple of elements in A^(x)i
            self.right = right or tuple()  # the tuple of elements in A^(x)j

        # converts this generator to an actual element
        def to_element(self) -> Bimodule.Element:
            return Bimodule.Element(self.module, len(self.left), len(self.right), {self: self.module.ring.one()})

        # returns m
        def get_module_generator(self):
            return Bimodule.Generator(self.module, self.key,
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
        def __add__(self, other) -> Bimodule.Element:
            return self.to_element() + other

        # tensor product
        @multimethod
        def __pow__(self, other: AMinus.Generator) -> Bimodule.Generator:
            assert self.rightmost_idempotent() == other.left_idempotent()
            return Bimodule.Generator(self.module, self.key, self.left_idempotent, self.right_idempotent,
                                      self.left, self.right + (other,))

        # tensor product
        @multimethod
        def __pow__(self, other: AMinus.Element) -> Bimodule.Element:
            return self.to_element() ** other

        @multimethod
        def __pow__(self, other) -> Bimodule.Generator:
            out = self
            for gen in other:
                out = out ** gen
            return out

        # tensor product
        @multimethod
        def __rpow__(self, other: AMinus.Generator) -> Bimodule.Generator:
            assert other.right_idempotent() == self.leftmost_idempotent()
            return Bimodule.Generator(self.module, self.key, self.left_idempotent, self.right_idempotent,
                                      (other,) + self.left, self.right)

        @multimethod
        def __rpow__(self, other: AMinus.Element) -> Bimodule.Element:
            return other ** self.to_element()

        @multimethod
        def __rpow__(self, other) -> Bimodule.Generator:
            out = self
            for gen in reversed(other):
                out = gen ** out
            return out

        # scalar multiplication
        def __rmul__(self, other: Z2Polynomial) -> Bimodule.Element:
            return other * self.to_element()

        def __str__(self):
            return str((self.left, self.key, self.right))

        def __repr__(self):
            return str((self.left, self.key, self.right))

        @multimethod
        def __eq__(self, other: Bimodule.Generator):
            return self.module == other.module and \
                   self.key == other.key and \
                   self.left == other.left and \
                   self.right == other.right

        @multimethod
        def __eq__(self, other):
            return self.to_element() == other

        def __hash__(self):
            return hash((self.module, self.key, self.left, self.right))


# represents a type DA bimodule
class TypeDA(Bimodule):
    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 left_scalar_action: Z2PolynomialRing.Map, right_scalar_action: Z2PolynomialRing.Map,
                 graph: MultiDiGraph=None):
        super().__init__(ring, left_algebra, right_algebra, left_scalar_action, right_scalar_action, graph=graph)

    # add the given generator to this module
    def add_generator(self, generator: Bimodule.Generator) -> None:
        self.graph.add_node(generator)

    # add the structure map (input |-> output) to this module
    def add_structure_map(self, input: Bimodule.Generator, output: Bimodule.Element) -> None:
        assert len(input.left) == output.j == 0 and output.i == 1
        x = input.get_module_generator()
        right_gens = input.right
        for gen_out, c_out in output.coefficients.items():
            y = gen_out.get_module_generator()
            for left_gen, c in gen_out.left[0].coefficients.items():
                current = self.graph.get_edge_data(x, y, key=(left_gen, right_gens))
                if current is None:
                    self.graph.add_edge(x, y, key=(left_gen, right_gens), c=self.ring.zero())
                    current = self.graph.get_edge_data(x, y, key=(left_gen, right_gens))
                current['c'] += c_out * self.left_scalar_action.apply(c)

    # turns this bimodule into a graphviz-compatible format
    def to_agraph(self, idempotents=True) -> AGraph:
        graph = AGraph(strict=False, directed=True)
        for generator in self.graph.nodes:
            graph.add_node(generator.key,
                           shape='box',
                           fontname='Arial')
        for x, y, (left, right), d in self.graph.edges(keys=True, data=True):
            c = d['c']
            if not idempotents and len(right) == 1 and right[0].is_idempotent():
                continue
            graph.add_edge(x.key, y.key,
                           label=str((left, c, right)),
                           dir='forward',
                           color=['black', 'blue', 'red', 'green'][len(right)],
                           fontname='Arial')
        graph.layout('dot')
        return graph

    # reduce this bimodule using the cancellation lemma
    def reduced(self) -> TypeDA:
        out = self
        reducible_edge = out.reducible_edge()
        while reducible_edge is not None:
            out = out.reduce_edge(*reducible_edge)
            reducible_edge = out.reducible_edge()
        return out

    def reducible_edge(self) -> Optional[Tuple]:
        for x in self.graph:
            for y in self.graph[x]:
                if len(self.graph[x][y]) == 1:
                    k, d = list(self.graph[x][y].items())[0]
                    if k[0].is_idempotent() and k[1] == tuple() and d['c'] == self.ring.one():
                        return x, y, k, d
        return None

    def reduce_edge(self, x, y, k, d) -> TypeDA:
        reduced_graph = MultiDiGraph(
            self.graph.subgraph([node for node in self.graph.nodes if node != x and node != y]))
        reduced_module = TypeDA(self.ring, self.left_algebra, self.right_algebra,
                                self.left_scalar_action, self.right_scalar_action, reduced_graph)

        left = k[0]
        right = k[1]
        c = d['c']

        for w, _, (left_wy, right_wy), d_wy in self.graph.in_edges(y, keys=True, data=True):
            if w == x or w == y:
                continue
            c_wy = d_wy['c']
            for _, z, (left_xz, right_xz), d_xz in self.graph.out_edges(x, keys=True, data=True):
                if z == x or z == y:
                    continue
                c_xz = d_xz['c']
                reduced_module.add_structure_map(w ** (right_wy + right_xz),
                                                 c_wy * c_xz * ((left_wy * left * left_xz) ** z))

        return reduced_module

    # tensor product of type DA structures
    # assumes self is bounded, other may or may not be
    def __pow__(self, other: TypeDA) -> TypeDA:
        assert self.right_algebra.ss == other.left_algebra.ss

        in_m, in_n = self.ring.tensor_inclusions(other.ring)

        out = TypeDA(in_m.target, self.left_algebra, other.right_algebra,
                     in_m.compose(self.left_scalar_action), in_n.compose(other.right_scalar_action))

        for x_m in self.graph.nodes:
            for x_n in other.graph.nodes:
                if x_m.right_idempotent == x_n.left_idempotent:
                    out.add_generator(Bimodule.Generator(out, (x_m.key, x_n.key),
                                                         x_m.left_idempotent, x_n.right_idempotent))

        for x_m in self.graph.nodes:
            for x_n in other.graph.nodes:
                if x_m.right_idempotent != x_n.left_idempotent:
                    continue
                x = Bimodule.Generator(out, (x_m.key, x_n.key), x_m.left_idempotent, x_n.right_idempotent)
                for _, y_m, (left_m, right_m), d_m in self.graph.out_edges(x_m, keys=True, data=True):
                    for y_n in other.graph.nodes:
                        if y_m.right_idempotent != y_n.left_idempotent:
                            continue
                        c_m = d_m['c']
                        y = Bimodule.Generator(out, (y_m.key, y_n.key), y_m.left_idempotent, y_n.right_idempotent)
                        for left_n, right_n, c_n in other.delta_n(len(right_m), x_n, y_n):
                            if left_n != right_m:
                                continue
                            out.add_structure_map(x ** right_n, left_m ** (in_m.apply(c_m) * in_n.apply(c_n) * y))

        return out

    # returns [left, right, coefficient)] representing the delta_n paths from source to target
    def delta_n(self, n, source, target) -> List[Tuple, Tuple, Z2Polynomial]:
        if n == 0:
            return [(tuple(), tuple(), self.ring.one())]
        else:
            out = []
            for new_target, _, k, d in self.graph.in_edges(target, keys=True, data=True):
                left = k[0]
                right = k[1]
                c = d['c']
                out += [(more_left + (left,), right, more_c * c)
                        for more_left, more_c in self.delta_n_helper(n - 1, source, new_target, right)]
            return out

    # returns [(left, coefficient)] representing the delta_n paths from source ** right to target
    @lru_cache(maxsize=None)
    def delta_n_helper(self, n, source, target, current_right) -> List[Tuple, Z2Polynomial]:
        if n == 0:
            return [(tuple(), self.ring.one())]
        else:
            out = []
            for new_target, _, k, d in self.graph.in_edges(target, keys=True, data=True):
                left = k[0]
                right = k[1]
                c = d['c']
                if right != current_right:
                    continue
                out += [(more_left + (left,), more_c * c)
                        for more_left, more_c in self.delta_n_helper(n - 1, source, new_target, right)]
            return out
