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
class ChainComplex(Module):
    def __init__(self, ring: Z2PolynomialRing, graph: MultiDiGraph = None, gradings: dict = None):
        super().__init__(ring, None, None, None, None, graph, gradings)

    # add the structure map (input |-> output) to this module
    def add_structure_map(self, input: Module.TensorGenerator, output: Module.TensorElement) -> None:
        assert len(input.left) == len(input.right) == output.i == output.j == 0
        x = input
        for gen_out, c_out in output.coefficients.items():
            y = gen_out
            self.add_edge(x, y, (), c_out)

    def edge_is_reducible(self, x, y) -> bool:
        if x in self.graph and y in self.graph[x] and len(self.graph[x][y]) == 1:
            k, d = list(self.graph[x][y].items())[0]
            if d['c'] == self.ring.one():
                return True
        return False

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

    # turns this bimodule into a graphviz-compatible format
    def to_agraph(self) -> AGraph:
        graph = AGraph(strict=False, directed=True)
        for generator in self.graph.nodes:
            graph.add_node(str(generator.key)+str(self.gradings[generator]),
                           shape='box',
                           fontname='Arial')
        for x, y, (), d in self.graph.edges(keys=True, data=True):
            c = d['c']
            graph.add_edge(str(x.key)+str(self.gradings[x]), str(y.key)+str(self.gradings[y]),
                           label=str(c),
                           dir='forward',
                           color='black',
                           fontname='Arial')
        graph.layout('dot')
        return graph

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

    # returns the direct sum decomposition of this module
    def decomposed(self) -> List[ChainComplex]:
        return [ChainComplex(self.ring, MultiDiGraph(self.graph.subgraph(component)),
                             {g: self.gradings[g] for g in self.graph.subgraph(component).nodes})
                for component in nx.weakly_connected_components(self.graph)]

    @staticmethod
    def direct_sum(modules: List) -> ChainComplex:
        new_graph = nx.union_all([da.graph for da in modules])
        return ChainComplex(modules[0].ring, new_graph,
                            {g: da.gradings[g] for da in modules for g in da.gradings.keys()})

    def reduce_edge(self, x, y, k, d) -> None:
        assert self.edge_is_reducible(x, y)

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
