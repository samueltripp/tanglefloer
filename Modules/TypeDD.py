from __future__ import annotations

from networkx import MultiDiGraph

from Modules.Module import Module
from SignAlgebra.AMinus import AMinus
from SignAlgebra.Z2PolynomialRing import Z2PolynomialRing


# represents a type DD bimodule
class TypeDD(Module):
    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 left_scalar_action: Z2PolynomialRing.Map, right_scalar_action: Z2PolynomialRing.Map,
                 graph: MultiDiGraph = None):
        super().__init__(ring, left_algebra, right_algebra, left_scalar_action, right_scalar_action, graph)

    @staticmethod
    def valid_input_gen(g):
        return g.num_left_factors() == 0 and g.num_right_factors() == 0

    @staticmethod
    def valid_output_gen(g):
        return g.num_left_factors() == 1 and g.num_right_factors() == 1

    @staticmethod
    def is_idempotent_edge_data(left, c, right):
        return False

    @staticmethod
    def edge_color(left, c, right):
        return 'black'

    @staticmethod
    def edge_is_reducible(left, c, right) -> bool:
        return left.to_algebra().is_idempotent() and right.to_algebra().is_idempotent() \
               and c == c.ring.one()

    def reduce_edge(self, x, y, k, d) -> None:
        in_edges = list(self.graph.in_edges(y, keys=True, data=True))
        out_edges = list(self.graph.out_edges(x, keys=True, data=True))

        self.graph.remove_nodes_from([x, y])

        left = k[0]
        right = k[1]

        for w, _, (left_wy, right_wy), d_wy in in_edges:
            if w == x or w == y:
                continue
            c_wy = d_wy['c']
            for _, z, (left_xz, right_xz), d_xz in out_edges:
                if z == x or z == y:
                    continue
                c_xz = d_xz['c']
                self.add_structure_map(w,
                                       (left_wy.to_algebra() * left.to_algebra() * left_xz.to_algebra())
                                       ** (c_wy * c_xz * z)
                                       ** (right_wy.to_algebra() * right.to_algebra() * right_xz.to_algebra()))
