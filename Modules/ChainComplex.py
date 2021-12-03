from __future__ import annotations

from typing import List

from multimethod import multimethod
from networkx import MultiDiGraph

from Modules.Module import Module
from SignAlgebra.AMinus import AMinus
from SignAlgebra.Z2PolynomialRing import Z2PolynomialRing


class ChainComplex(Module):

    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 left_scalar_action: Z2PolynomialRing.Map, right_scalar_action: Z2PolynomialRing.Map,
                 graph: MultiDiGraph = None):
        super().__init__(ring, left_algebra, right_algebra, left_scalar_action, right_scalar_action, graph)

    @staticmethod
    def valid_input_gen(g):
        return g.num_left_factors() == 0 and g.num_right_factors() == 0

    @staticmethod
    def valid_output_gen(g):
        return g.num_left_factors() == 0 and g.num_right_factors() == 0

    @staticmethod
    def is_idempotent_edge_data(left, c, right):
        return False

    @staticmethod
    def edge_color(left, c, right):
        return 'black'

    @multimethod
    def d(self, elt: Module.TensorElement) -> Module.TensorElement:
        out = self.zero()

        for g, c in elt.coefficients.items():
            out += c * self.d(g)

        return out

    @multimethod
    def d(self, x: Module.TensorGenerator) -> Module.TensorElement:
        out = self.zero()

        for _, y, (), d in self.graph.out_edges(x, keys=True, data=True):
            out += d['c'] * y

        return out

    def d_squared_is_zero(self) -> bool:
        for x in self.graph.nodes:
            if self.d(self.d(x)) != self.zero():
                print(x)
                print(self.d(x))
                print(self.d(self.d(x)))
                return False
        return True

    def m2_def(self) -> List[str]:
        arrows_per_def = 50

        gens = list(self.graph.nodes)

        arrow_strings = []
        for i, x in enumerate(gens):
            for j, y in enumerate(gens):
                if x == y:
                    continue
                if y in self.graph[x]:
                    k, d = list(self.graph[x][y].items())[0]
                    c = d['c']
                    arrow_strings += [f"({j},{i}) => {c}"]

        out = [f"R = ZZ/2[{','.join(self.ring.variables)}]"]
        out += [f"M = R^{len(gens)}"]
        out += [f"d = map(M, M, 0)"]
        for i in range(0, (len(arrow_strings) // arrows_per_def) + 1):
            out += [f"d = d + map(M, M, {{{', '.join(arrow_strings[i * arrows_per_def:(i+1) * arrows_per_def])}}})"]
        out += ["trim' = Q -> image(generators Q)/intersect(image generators Q, image relations Q)"]
        out += ["H = trim' trim ((ker d) / (image d))"]

        return out

    def write_m2_def(self, filename: str):
        with open(filename, 'w') as out:
            out.writelines([line+'\n' for line in self.m2_def()])

    @staticmethod
    def edge_is_reducible(left, c, right) -> bool:
        return left == left.tensor_algebra.one() \
               and right == right.tensor_algebra.one() \
               and c == c.ring.one()

    def reduce_edge(self, x, y, k, d) -> None:
        in_edges = list(self.graph.in_edges(y, keys=True, data=True))
        out_edges = list(self.graph.out_edges(x, keys=True, data=True))

        self.graph.remove_nodes_from([x, y])

        for w, _, (), d_wy in in_edges:
            if w == x or w == y:
                continue
            c_wy = d_wy['c']
            for _, z, (), d_xz in out_edges:
                if z == x or z == y:
                    continue
                c_xz = d_xz['c']
                self.add_structure_map(w, c_wy * c_xz * z)
