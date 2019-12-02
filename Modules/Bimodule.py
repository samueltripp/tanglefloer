from __future__ import annotations
from networkx import MultiDiGraph
from typing import Iterable, Any

from pygraphviz import AGraph

from Modules import ETangleStrands
from SignAlgebra.AMinus import AMinus
from Modules.CTMinus import *


# Base class for Type DD, AA, DA, and AD structures
class Bimodule:
    def __init__(self, left_algebra: AMinus, right_algebra: AMinus,
                 generators: Iterable[Bimodule.Generator], maps: Iterable[Bimodule.Edge]):
        self.left_algebra = left_algebra
        self.right_algebra = right_algebra

        self.graph = MultiDiGraph()
        for gen in generators:
            self.graph.add_node(gen.key, left_idempotent=gen.left_idempotent, right_idempotent=gen.right_idempotent)
        for edge in maps:
            if edge.c != edge.target_diagram.etangle.polyring.zero():
                self.graph.add_edge(edge.source_diagram, edge.target_diagram,
                                    c=edge.c, left=edge.left, right=edge.right)

    @classmethod
    def from_strand_diagrams(cls, left_algebra: AMinus, right_algebra: AMinus,
                             generators: Iterable[ETangleStrands], maps: Iterable[Bimodule.Edge]):
        return Bimodule(left_algebra, right_algebra,
                        [Bimodule.Generator(sd, sd.left_idempotent(), sd.right_idempotent())
                         for sd in generators], maps)

    def __repr__(self) -> str:
        return str(self.__dict__)

    def to_agraph(self):
        out = AGraph(size='10,10', nodesep='1')
        for node in self.graph.nodes:
            out.add_node(node, label=str(dict(node.left_strands)) + str(dict(node.right_strands)), shape='box')
        for source in self.graph:
            for target in self.graph[source]:
                for i in self.graph[source][target]:
                    edge = self.graph[source][target][i]
                    out.add_edge(source, target, label=str((edge['left'], edge['c'], edge['right'])), dir='forward')
        out.layout('dot')
        return out

    class Generator:
        def __init__(self, key, left_idempotent: AMinus.Element, right_idempotent: AMinus.Element):
            self.key = key
            self.left_idempotent = left_idempotent
            self.right_idempotent = right_idempotent

    class Element:
        # d - {StrandDiagram: Z2Polynomial}
        def __init__(self, d=None):
            if d is None:
                d = {}
            self.d = {}
            for sd, c in d.items():
                if c != sd.etangle.polyring.zero():
                    self.d[sd] = c

        def __add__(self, other: Bimodule.Element) -> Bimodule.Element:
            out_d = dict(self.d)
            for sd in other.d:
                if sd in self.d:
                    out_d[sd] = self.d[sd] + other.d[sd]
                else:
                    out_d[sd] = other.d[sd]
            return Bimodule.Element(out_d)

        def __rmul__(self, other: Z2Polynomial):
            d_out = dict(self.d)
            for k in d_out:
                d_out[k] = other * d_out[k]
            return Bimodule.Element(d_out)

        def __eq__(self, other: Bimodule.Element) -> bool:
            return self.d == other.d

        def __repr__(self) -> str:
            return str(self.d)

    class Edge:
        def __init__(self, source_diagram, target_diagram, c: Z2Polynomial, left: Tuple, right: Tuple):
            self.source_diagram = source_diagram
            self.target_diagram = target_diagram
            self.c = c
            self.left = left
            self.right = right

        def __repr__(self) -> str:
            return str(self.__dict__)

        def __eq__(self, other: Bimodule.Edge) -> bool:
            return self.source_diagram == other.source_diagram and \
                   self.target_diagram == other.target_diagram and \
                   self.c == other.c and \
                   self.left == other.left and \
                   self.right == other.right


class TypeDA(Bimodule):
    def __init__(self, left_algebra: AMinus, right_algebra: AMinus,
                 generators: Iterable[Bimodule.Generator], maps: Iterable[Bimodule.Edge]):
        super().__init__(left_algebra, right_algebra, generators, maps)

    # TODO: this won't work because of polynomial ring mismatch
    def tensor(self, other: TypeDA) -> TypeDA:
        assert self.right_algebra == other.left_algebra

        generators = [Bimodule.Generator((xm, xn), xm['left_idempotent'], xn['right_idempotent'])
                      for xm in self.graph.nodes for xn in other.graph.nodes
                      if xm['right_idempotent'] == xn['left_idempotent']]

        maps = []
        for (xm, xn) in generators:
            for ym in self.graph[xm]:
                for i in self.graph[xm][ym]:
                    delta_1 = self.graph[xm][ym][i]
                    if len(delta_1['left']) > 1:
                        continue
                    for yn in self.graph[xn]:
                        for j in self.graph[xn][yn]:
                            delta_n = self.graph[xn][yn][j]
                            if delta_1['right'] != delta_n['left']:
                                continue
                            maps.append(
                                Bimodule.Edge((xm, xn), (ym, yn),
                                              delta_1['c'] * delta_n['c'], delta_1['left'], delta_n['right']))

        return TypeDA(self.left_algebra, other.right_algebra, generators, maps)
