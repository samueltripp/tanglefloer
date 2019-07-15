from __future__ import annotations
from networkx import MultiDiGraph
from typing import Iterable, Any
from SignAlgebra.AMinus import *
import itertools as it


# Base class for Type DD, AA, DA, and AD structures
class Bimodule:
    def __init__(self, left_algebra: AMinus, right_algebra: AMinus,
                 generators: Iterable, maps: Iterable[Bimodule.Edge]):
        self.left_algebra = left_algebra
        self.right_algebra = right_algebra

        self.graph = MultiDiGraph()
        self.graph.add_nodes_from(generators)
        for edge in maps:
            if edge.c != 0:
                self.graph.add_edge(edge.source, edge.target, c=edge.c, left=edge.left, right=edge.right)

    # I wrote this for DD tensor AA but it might be generic enough for all cases
    # the only specialized parts are the variable names delta and m
    def tensor(self, other: Bimodule) -> Bimodule:
        assert self.right_algebra == other.left_algebra

        generators = it.product(self.graph.nodes, other.graph.nodes)

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

    class Edge:
        def __init__(self, source, target, c: Z2Polynomial, left: AMinusElement, right: AMinusElement):
            self.source = source
            self.target = target
            self.c = c
            self.left = left
            self.right = right
