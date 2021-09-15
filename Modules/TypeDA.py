from __future__ import annotations

from collections import defaultdict
from functools import lru_cache

from networkx import MultiDiGraph
import networkx as nx
from typing import Iterable
from pygraphviz import AGraph
from Modules import ETangleStrands
from Modules.ChainComplex import ChainComplex
from SignAlgebra.AMinus import AMinus
from Modules.CTMinus import *
from multimethod import *
from frozendict import *
from SignAlgebra.Z2PolynomialRing import *

from Modules.Module import Module


# represents a type DA bimodule
class TypeDA(Module):
    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 left_scalar_action: Z2PolynomialRing.Map, right_scalar_action: Z2PolynomialRing.Map,
                 graph: MultiDiGraph = None):
        super().__init__(ring, left_algebra, right_algebra, left_scalar_action, right_scalar_action, graph)

    # add the structure map (input |-> output) to this module
    def add_structure_map(self, input: Module.TensorGenerator, output: Module.TensorElement) -> None:
        assert len(input.left) == output.j == 0 and output.i == 1
        x = input.get_module_generator()
        right_gens = input.right
        for gen_out, c_out in output.coefficients.items():
            y = gen_out.get_module_generator()
            left_gen = gen_out.left[0]
            self.add_edge(x, y, (left_gen, right_gens), c_out)

    def edge_is_reducible(self, x, y) -> bool:
        if x in self.graph and y in self.graph[x] and len(self.graph[x][y]) == 1:
            k, d = list(self.graph[x][y].items())[0]
            if k[0].is_idempotent() and k[1] == tuple() and d['c'] == self.ring.one():
                return True
        return False

    # turns this bimodule into a graphviz-compatible format
    def to_agraph(self, idempotents=True) -> AGraph:
        graph = AGraph(strict=False, directed=True)
        for generator, data in self.graph.nodes(data=True):
            graph.add_node(str(generator.key) + '[' + str(data['grading'])[1:-1] + ']',
                           shape='box',
                           fontname='Arial')
        for x, y, (left, right), d in self.graph.edges(keys=True, data=True):
            c = d['c']
            x_grading = self.graph.nodes(data=True)[x]['grading']
            y_grading = self.graph.nodes(data=True)[y]['grading']
            if not idempotents and len(right) == 1 and right[0].is_idempotent():
                continue
            graph.add_edge(str(x.key) + '[' + str(x_grading)[1:-1] + ']',
                           str(y.key) + '[' + str(y_grading)[1:-1] + ']',
                           label=' '+str((left, c, right))+' ',
                           dir='forward',
                           color=['black', 'blue', 'red', 'green', 'purple'][min(len(right), 4)],
                           fontname='Arial',
                           decorate=True)
        graph.layout('dot')
        return graph

    # probably wrong, don't actually know if A^- can be augmented
    def to_chain_complex(self) -> ChainComplex:
        out = ChainComplex(self.ring)

        for generator, generator_data in self.graph.nodes(data=True):
            out.add_generator(generator, generator_data['grading'])

        for x, y, (left, right), d in self.graph.edges(keys=True, data=True):
            if left.is_idempotent() and len(right) == 0:
                out.add_structure_map(x, d['c'] * y)

        return out

    def reduce_edge(self, x, y, k, d) -> None:
        assert self.edge_is_reducible(x, y)

        in_edges = list(self.graph.in_edges(y, keys=True, data=True))
        out_edges = list(self.graph.out_edges(x, keys=True, data=True))

        self.graph.remove_nodes_from([x, y])

        left = k[0]

        for w, _, (left_wy, right_wy), d_wy in in_edges:
            if w == x or w == y:
                continue
            c_wy = d_wy['c']
            for _, z, (left_xz, right_xz), d_xz in out_edges:
                if z == x or z == y:
                    continue
                c_xz = d_xz['c']
                self.add_structure_map(
                    w ** (right_wy + right_xz),
                    c_wy * c_xz * ((left_wy * left * left_xz) ** z))

    # tensor product of type DA structures
    # assumes self is bounded, other may or may not be
    def __pow__(self, other: TypeDA) -> TypeDA:
        assert self.right_algebra == other.left_algebra

        left_right_iso = Z2PolynomialRing.Map.identity(other.left_algebra.ring, self.right_algebra.ring)
        in_m, in_n = self.right_scalar_action.compose(left_right_iso).pushout_inclusions(other.left_scalar_action)

        out = TypeDA(in_m.target, self.left_algebra, other.right_algebra,
                     in_m.compose(self.left_scalar_action), in_n.compose(other.right_scalar_action))

        for x_m, x_m_data in self.graph.nodes(data=True):
            for x_n, x_n_data in other.graph.nodes(data=True):
                x_m_grading = x_m_data['grading']
                x_n_grading = x_n_data['grading']
                if x_m.right_idempotent == x_n.left_idempotent:
                    out.add_generator(Module.TensorGenerator(out, (x_m.key, x_n.key),
                                                             x_m.left_idempotent, x_n.right_idempotent),
                                      (x_m_grading[0] + x_n_grading[0], x_m_grading[1] + x_n_grading[1]))

        for x_m in self.graph.nodes:
            for x_n in other.graph.nodes:
                if x_m.right_idempotent != x_n.left_idempotent:
                    continue
                x = Module.TensorGenerator(out, (x_m.key, x_n.key),
                                           x_m.left_idempotent, x_n.right_idempotent)
                for _, y_m, (left_m, right_m), d_m in self.graph.out_edges(x_m, keys=True, data=True):
                    for right_n, y_n, c_n in other.delta_n(right_m, x_n):
                        if y_m.right_idempotent != y_n.left_idempotent:
                            continue
                        c_m = d_m['c']
                        y = Module.TensorGenerator(out, (y_m.key, y_n.key),
                                                   y_m.left_idempotent, y_n.right_idempotent)
                        out.add_structure_map(
                            x ** right_n,
                            left_m ** (in_m.apply(c_m) * in_n.apply(c_n) * y))

        # assert (out.to_chain_complex().d_squared_is_zero())  # probably slow, take out when not debugging

        return out

    # returns [(right, target, coefficient)]
    #   representing the delta_n paths starting at source outputting left
    def delta_n(self, left, source) -> List[Tuple[Tuple, Module.TensorGenerator, Z2Polynomial]]:
        if len(left) == 0:
            return [(tuple(), source, self.ring.one())]
        else:
            out = []
            for _, new_source, k, d in self.graph.out_edges(source, keys=True, data=True):
                if k[0] != left[0]:
                    continue
                right = k[1]
                c = d['c']
                out += [(right + more_right, target, more_c * c)
                        for more_right, target, more_c in
                        self.delta_n(left[1:], new_source)]
            return out
