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

from Modules.Module import Module


# represents a type DA bimodule
class TypeDA(Module):
    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 right_scalar_action: Z2PolynomialRing.Map,
                 graph: MultiDiGraph = None):
        super().__init__(ring, left_algebra, right_algebra, None, right_scalar_action, graph=graph)

    # add the structure map (input |-> output) to this module
    def add_structure_map(self, input: Module.TensorGenerator, output: Module.TensorElement) -> None:
        assert len(input.left) == output.j == 0 and output.i == 1
        x = input.get_module_generator()
        right_gens = input.right
        for gen_out, c_out in output.coefficients.items():
            y = gen_out.get_module_generator()
            left_gen = gen_out.left[0]
            left_monomial = gen_out.left_monomial
            current = self.graph.get_edge_data(x, y, key=(left_monomial, left_gen, right_gens))
            if current is None:
                self.graph.add_edge(x, y, key=(left_monomial, left_gen, right_gens), c=self.ring.zero())
                current = self.graph.get_edge_data(x, y, key=(left_monomial, left_gen, right_gens))
            current['c'] += c_out
            if current['c'] == self.ring.zero():
                self.graph.remove_edge(x, y, key=(left_monomial, left_gen, right_gens))

    # turns this bimodule into a graphviz-compatible format
    def to_agraph(self, idempotents=True) -> AGraph:
        graph = AGraph(strict=False, directed=True)
        for generator in self.graph.nodes:
            graph.add_node(generator.key,
                           shape='box',
                           fontname='Arial')
        for x, y, (left_monomial, left, right), d in self.graph.edges(keys=True, data=True):
            c = d['c']
            if not idempotents and len(right) == 1 and right[0].is_idempotent():
                continue
            graph.add_edge(x.key, y.key,
                           label=str((left_monomial, left, c, right)),
                           dir='forward',
                           color=['black', 'blue', 'red', 'green', 'purple'][len(right)],
                           fontname='Arial')
        graph.layout('dot')
        return graph

    # returns the direct sum decomposition of this module
    def decomposed(self) -> List[TypeDA]:
        return [TypeDA(self.ring, self.left_algebra, self.right_algebra,
                       self.right_scalar_action, MultiDiGraph(self.graph.subgraph(component)))
                for component in nx.weakly_connected_components(self.graph)]

    @staticmethod
    def direct_sum(modules: List) -> TypeDA:
        new_graph = nx.union_all([da.graph for da in modules])
        return TypeDA(modules[0].ring, modules[0].left_algebra, modules[0].right_algebra,
                      modules[0].right_scalar_action, new_graph)

    # O(V^2)
    def reducible_edge(self) -> Optional[Tuple]:
        for x in self.graph:
            for y in self.graph[x]:
                if len(self.graph[x][y]) == 1:
                    k, d = list(self.graph[x][y].items())[0]
                    if k[0].to_polynomial() == self.left_algebra.ring.one() and \
                            k[1].is_idempotent() and k[2] == tuple() and d['c'] == self.ring.one():
                        return x, y, k, d
        return None

    def reduce_edge(self, x, y, k, d) -> None:
        in_edges = list(self.graph.in_edges(y, keys=True, data=True))
        out_edges = list(self.graph.out_edges(x, keys=True, data=True))

        self.graph.remove_nodes_from([x, y])

        left = k[1]

        for w, _, (left_monomial_wy, left_wy, right_wy), d_wy in in_edges:
            if w == x or w == y:
                continue
            c_wy = d_wy['c']
            for _, z, (left_monomial_xz, left_xz, right_xz), d_xz in out_edges:
                if z == x or z == y:
                    continue
                c_xz = d_xz['c']
                self.add_structure_map(
                    w ** (right_wy + right_xz),
                    c_wy * c_xz *
                    ((left_monomial_wy.to_polynomial() * left_monomial_xz.to_polynomial() * left_wy * left * left_xz)
                     ** z))

    # tensor product of type DA structures
    # assumes self is bounded, other may or may not be
    def __pow__(self, other: TypeDA) -> TypeDA:
        assert self.right_algebra.ss == other.left_algebra.ss

        in_m, in_n = self.ring.tensor_inclusions(other.ring)
        scalar_map = Z2PolynomialRing.Map.identity(other.left_algebra.ring, self.right_algebra.ring)

        out = TypeDA(in_m.target, self.left_algebra, other.right_algebra, in_n.compose(other.right_scalar_action))

        for x_m in self.graph.nodes:
            for x_n in other.graph.nodes:
                if x_m.right_idempotent == x_n.left_idempotent:
                    out.add_generator(Module.TensorGenerator(out, (x_m.key, x_n.key),
                                                             x_m.left_idempotent, x_n.right_idempotent))

        for x_m in self.graph.nodes:
            for x_n in other.graph.nodes:
                if x_m.right_idempotent != x_n.left_idempotent:
                    continue
                x = Module.TensorGenerator(out, (x_m.key, x_n.key),
                                           x_m.left_idempotent, x_n.right_idempotent)
                for _, y_m, (left_monomial, left_m, right_m), d_m in self.graph.out_edges(x_m, keys=True, data=True):
                    for y_n in other.graph.nodes:
                        if y_m.right_idempotent != y_n.left_idempotent:
                            continue
                        c_m = d_m['c']
                        y = Module.TensorGenerator(out, (y_m.key, y_n.key),
                                                   y_m.left_idempotent, y_n.right_idempotent)
                        for left_monomial_n, left_n, right_n, c_n in other.delta_n(len(right_m), x_n, y_n):
                            if left_n != right_m:
                                continue
                            out.add_structure_map(
                                x ** right_n,
                                (left_monomial.to_polynomial() * left_m) **
                                (in_m.apply(c_m * self.right_scalar_action.apply(
                                    scalar_map.apply(left_monomial_n.to_polynomial()))) *
                                 in_n.apply(c_n) * y))

        return out

    # returns [(left_monomial, left, right, coefficient)] representing the delta_n paths from source to target
    def delta_n(self, n, source, target) -> List[Tuple, Tuple, Z2Polynomial]:
        if n == 0:
            if source == target:
                return [(Z2Monomial(self.left_algebra.ring, {}), tuple(), tuple(), self.ring.one())]
            else:
                return []
        else:
            out = []
            for new_target, _, k, d in self.graph.in_edges(target, keys=True, data=True):
                left_monomial = k[0]
                left = k[1]
                right = k[2]
                c = d['c']
                out += [(more_left_monomial * left_monomial, more_left + (left,), right, more_c * c)
                        for more_left_monomial, more_left, more_c in
                        self.delta_n_helper(n - 1, source, new_target, right)]
            return out

    # returns [(left_monomial, left, coefficient)] representing the delta_n paths from source ** right to target
    @lru_cache(maxsize=None)
    def delta_n_helper(self, n, source, target, current_right) -> List[Tuple, Z2Polynomial]:
        if n == 0:
            if source == target:
                return [(Z2Monomial(self.left_algebra.ring, {}), tuple(), self.ring.one())]
            else:
                return []
        else:
            out = []
            for new_target, _, k, d in self.graph.in_edges(target, keys=True, data=True):
                left_monomial = k[0]
                left = k[1]
                right = k[2]
                c = d['c']
                if right != current_right:
                    continue
                out += [(more_left_monomial * left_monomial, more_left + (left,), more_c * c)
                        for more_left_monomial, more_left, more_c
                        in self.delta_n_helper(n - 1, source, new_target, right)]
            return out
