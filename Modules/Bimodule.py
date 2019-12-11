from __future__ import annotations
from networkx import MultiDiGraph
from typing import Iterable, Any
from pygraphviz import AGraph
from Modules import ETangleStrands
from SignAlgebra.AMinus import AMinus
from Modules.CTMinus import *


# Base class for Type DD, AA, DA, and AD structures
class Bimodule:
    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 generators: Iterable[Bimodule.Generator], maps: Iterable[Bimodule.Edge]):
        self.ring = ring
        self.left_algebra = left_algebra
        self.right_algebra = right_algebra

        self.graph = MultiDiGraph()
        for gen in generators:
            self.graph.add_node(gen)
        for edge in maps:
            if edge.c != ring.zero():
                self.graph.add_edge(edge.source, edge.target,
                                    c=edge.c, left=edge.left, right=edge.right)

    def __repr__(self) -> str:
        return str(self.__dict__)

    def to_agraph(self, idempotents=True):
        out = AGraph(strict=False, directed=True)
        for node in self.graph.nodes:
            out.add_node(node,
                         shape='box',
                         fontname='Arial')
        for source in self.graph:
            for target in self.graph[source]:
                if not idempotents and target == source:
                    continue
                for i in self.graph[source][target]:
                    edge = self.graph[source][target][i]
                    color = 'black' if edge['right'] == tuple() else 'blue'
                    out.add_edge(source, target,
                                 label=str((edge['left'], edge['c'], edge['right'])),
                                 dir='forward',
                                 color=color,
                                 fontname='Arial')
        out.layout('dot')
        return out

    class Generator:
        def __init__(self, key, left_idempotent: AMinus.Element, right_idempotent: AMinus.Element):
            self.key = key
            self.left_idempotent = left_idempotent
            self.right_idempotent = right_idempotent

        def __eq__(self, other: Bimodule.Generator):
            return self.key == other.key and \
                   self.left_idempotent == other.left_idempotent and \
                   self.right_idempotent == other.right_idempotent

        def __str__(self):
            return str(self.key)

        def __repr__(self):
            return str(self.key)

        def __hash__(self):
            return hash(self.key)

    class Element:
        # d - {Bimodule.Generator: Z2Polynomial}
        def __init__(self, d=None):
            if d is None:
                d = {}
            self.d = {}
            for g, c in d.items():
                if c != c.ring.zero():
                    self.d[g] = c

        def __add__(self, other: Bimodule.Element) -> Bimodule.Element:
            out_d = dict(self.d)
            for g in other.d:
                if g in self.d:
                    out_d[g] = self.d[g] + other.d[g]
                else:
                    out_d[g] = other.d[g]
            return Bimodule.Element(out_d)

        def __rmul__(self, other: Z2Polynomial):
            d_out = dict(self.d)
            for g in d_out:
                d_out[g] = other * d_out[g]
            return Bimodule.Element(d_out)

        def __eq__(self, other: Bimodule.Element) -> bool:
            return self.d == other.d

        def __repr__(self) -> str:
            return str(self.d)

    class Edge:
        def __init__(self, source: Bimodule.Generator, target: Bimodule.Generator, c: Z2Polynomial, left: Tuple,
                     right: Tuple):
            self.source = source
            self.target = target
            self.c = c
            self.left = left
            self.right = right

        def __repr__(self) -> str:
            return str(self.__dict__)

        def __eq__(self, other: Bimodule.Edge) -> bool:
            return self.source == other.source and \
                   self.target == other.target and \
                   self.c == other.c and \
                   self.left == other.left and \
                   self.right == other.right

        def __hash__(self):
            return hash((self.source, self.target, self.c, self.left, self.right))


class TypeDA(Bimodule):
    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 generators: Iterable, maps: Iterable[Bimodule.Edge]):
        super().__init__(ring, left_algebra, right_algebra, generators, maps)

    def tensor(self, other: TypeDA) -> TypeDA:
        assert self.right_algebra.ss == other.left_algebra.ss

        generators = [Bimodule.Generator((xm, xn), xm.left_idempotent, xn.right_idempotent)
                      for xm in self.graph.nodes for xn in other.graph.nodes
                      if xm.right_idempotent == xn.left_idempotent]

        in1, in2 = self.ring.tensor_inclusions(other.ring)

        maps = set()
        for x in generators:
            xm = x.key[0]
            xn = x.key[1]
            for ym in self.graph[xm]:
                for i in self.graph[xm][ym]:
                    delta_1 = self.graph[xm][ym][i]
                    if len(delta_1['left']) > 1:
                        continue
                    for yn in other.graph[xn]:
                        if ym.right_idempotent != yn.left_idempotent:
                            continue
                        for j in other.graph[xn][yn]:
                            delta_n = other.graph[xn][yn][j]
                            if delta_1['right'] != delta_n['left']:
                                continue
                            maps.add(
                                Bimodule.Edge(
                                    Bimodule.Generator((xm, xn), xm.left_idempotent, xn.right_idempotent),
                                    Bimodule.Generator((ym, yn), ym.left_idempotent, yn.right_idempotent),
                                    in1.apply(delta_1['c']) * in2.apply(delta_n['c']), delta_1['left'], delta_n['right']))

        return TypeDA(in1.target, self.left_algebra, other.right_algebra, generators, maps)
