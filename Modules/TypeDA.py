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
from SignAlgebra.TensorAlgebra import *
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

    @staticmethod
    def valid_input_gen(g):
        return g.num_left_factors() == 0

    @staticmethod
    def valid_output_gen(g):
        return g.num_left_factors() == 1 and g.num_right_factors() == 0

    @staticmethod
    def is_idempotent_edge_data(left, c, right):
        return right.num_factors() == 1 and right.to_algebra().is_idempotent()

    @staticmethod
    def edge_color(left, c, right):
        return ['black', 'blue', 'red', 'green', 'purple'][min(right.num_factors(), 4)]

    @staticmethod
    def edge_is_reducible(left, c, right) -> bool:
        return left.to_algebra().is_idempotent() \
               and right == right.tensor_algebra.one() \
               and c == c.ring.one()

    def reduce_edge(self, x, y, k, d) -> None:
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
                    w ** (right_wy ** right_xz),
                    c_wy * c_xz * ((left_wy.to_algebra() * left.to_algebra() * left_xz.to_algebra()) ** z))

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

        return out

    # returns [(right, target, coefficient)]
    #   representing the delta_n paths starting at source outputting left
    def delta_n(self, left, source) -> List[Tuple[TensorAlgebra.Generator, Module.TensorGenerator, Z2Polynomial]]:
        if left.num_factors() == 0:
            return [(self.right_tensor_algebra.one_generator(), source, self.ring.one())]
        else:
            out = []
            for _, new_source, k, d in self.graph.out_edges(source, keys=True, data=True):
                if k[0].to_algebra() != left.factors[0]:
                    continue
                right = k[1]
                c = d['c']
                out += [(right ** more_right, target, more_c * c)
                        for more_right, target, more_c in
                        self.delta_n(TensorAlgebra.Generator(self.left_tensor_algebra, left.factors[1:]), new_source)]
            return out
