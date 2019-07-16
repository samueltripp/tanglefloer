from __future__ import annotations
from networkx import MultiDiGraph
from typing import Iterable, Any
from SignAlgebra.AMinus import *
from Modules.CTMinus import *
import itertools as it


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
            if edge.c != 0:
                self.graph.add_edge(edge.source_key, edge.target_key, c=edge.c, left=edge.left, right=edge.right)

    @classmethod
    def from_strand_diagrams(cls, left_algebra: AMinus, right_algebra: AMinus,
                             generators: Iterable[StrandDiagram], maps: Iterable[Bimodule.Edge]):
        return Bimodule(left_algebra, right_algebra,
                        [Bimodule.Generator(sd, sd.left_idempotent(), sd.right_idempotent())
                         for sd in generators], maps)

    # I wrote this for DD tensor AA but it might be generic enough for all cases
    # the only specialized parts are the variable names delta and m
    def tensor(self, other: Bimodule) -> Bimodule:
        assert self.right_algebra == other.left_algebra

        generators = [Bimodule.Generator((xm, xn), xm['left_idempotent'], xn['right_idempotent'])
                      for xm in self.graph.nodes for xn in other.graph.nodes
                      if xm['right_idempotent'] == xn['left_idempotent']]

        maps = set()
        for (xm, xn) in generators:
            for ym in self.graph[xm]:
                for i in self.graph[xm][ym]:
                    delta = self.graph[xm][ym][i]
                    for yn in self.graph[xn]:
                        for j in self.graph[xn][yn]:
                            m = self.graph[xn][yn][j]
                            if delta['right'] != m['left']:
                                continue
                            maps.add(Bimodule.Edge((xm, xn), (ym, yn), delta['c'] * m['c'], delta['left'], m['right']))

        return Bimodule(self.left_algebra, other.right_algebra, generators, maps)

    class Generator:
        def __init__(self, key, left_idempotent: AMinusElement, right_idempotent: AMinusElement):
            self.key = key
            self.left_idempotent = left_idempotent
            self.right_idempotent = right_idempotent

    class Edge:
        def __init__(self, source_key, target_key, c: Z2Polynomial, left: AMinusElement, right: AMinusElement):
            self.source_key = source_key
            self.target_key = target_key
            self.c = c
            self.left = left
            self.right = right
