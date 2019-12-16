from __future__ import annotations
from networkx import MultiDiGraph
from typing import Iterable, Any
from pygraphviz import AGraph
from Modules import ETangleStrands
from SignAlgebra.AMinus import AMinus
from Modules.CTMinus import *
from multimethod import *
from frozendict import *
from SignAlgebra.Z2PolynomialRing import *


# Base class for Type DD, AA, DA, and AD structures
class Bimodule:
    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 left_scalar_action: Z2PolynomialRing.Map, right_scalar_action: Z2PolynomialRing.Map):
        self.ring = ring
        self.left_algebra = left_algebra
        self.right_algebra = right_algebra
        self.left_scalar_action = left_scalar_action
        self.right_scalar_action = right_scalar_action
        self.generators = []
        self.structure_maps = {}  # {Bimodule.TensorGenerator: Bimodule.TensorElement}

    def __repr__(self) -> str:
        return str(self.__dict__)

    def add_generator(self, generator: Bimodule.TensorGenerator) -> None:
        self.generators.append(generator)

    def add_structure_map(self, input: Bimodule.TensorGenerator, output: Bimodule.TensorElement) -> None:
        if input in self.structure_maps:
            self.structure_maps[input] += output
        else:
            self.structure_maps[input] = output

    def generator(self, key, left_idempotent: AMinus.Element, right_idempotent: AMinus.Element) -> Bimodule.TensorGenerator:
        return Bimodule.TensorGenerator(self, key, tuple(), tuple())

    # returns the zero element of A^(x)i (x) M (x) A^(x)j
    def zero(self, i=0, j=0):
        return Bimodule.TensorElement(self, i, j, {})

    # represents an element of A^(x)i (x) M (x) A^(x)j
    class TensorElement:
        # coefficients - {TensorGenerator: Z2Polynomial}
        def __init__(self, module: Bimodule, i, j, coefficients=None):
            self.module = module
            self.i = i
            self.j = j
            if coefficients is None:
                coefficients = {}
            self.coefficients = frozendict(coefficients)
            self.simplify()

        def simplify(self):
            new_coefficients = dict(self.coefficients)
            for g, c in new_coefficients.items():
                if c == c.ring.zero():
                    del new_coefficients[g]
            self.coefficients = frozendict(new_coefficients)

        def __add__(self, other: Bimodule.TensorElement) -> Bimodule.TensorElement:
            assert self.i == other.i and self.j == other.j
            new_coefficients = dict(self.coefficients)
            for g in other.coefficients:
                if g in self.coefficients:
                    new_coefficients[g] = self.coefficients[g] + other.coefficients[g]
                else:
                    new_coefficients[g] = other.coefficients[g]
            return Bimodule.TensorElement(self.module, self.i, self.j, new_coefficients)

        @multimethod
        def __rmul__(self, other: Z2Polynomial) -> Bimodule.TensorElement:
            new_coefficients = dict(self.coefficients)
            for g in new_coefficients:
                new_coefficients[g] = other * new_coefficients[g]
            return Bimodule.TensorElement(self.module, self.i, self.j, new_coefficients)

        @multimethod
        def __rmul__(self, other: AMinus.Element) -> Bimodule.TensorElement:
            out = self.module.zero(self.i + 1, self.j)

            for g1, c1 in other.coefficients.items():
                for g2, c2 in self.coefficients.items():
                    out += (self.module.left_scalar_action.apply(c1) * c2) * (g1 * g2)

            return out

        def __mul__(self, other: AMinus.Element) -> Bimodule.TensorElement:
            out = self.module.zero(self.i, self.j + 1)

            for g1, c1 in self.coefficients.items():
                for g2, c2 in other.coefficients.items():
                    out += (c1 * self.module.right_scalar_action.apply(c2)) * (g1 * g2)

            return out

        def __eq__(self, other: Bimodule.TensorElement) -> bool:
            return self.coefficients == other.coefficients

        def __repr__(self) -> str:
            return str(self.coefficients)

    # represents a generator of A^(x)i (x) M (x) A^(x)j
    class TensorGenerator(TensorElement):
        def __init__(self, module: Bimodule, key,
                     left_elements: Tuple[AMinus.Generator, ...], right_elements: Tuple[AMinus.Generator, ...]):
            self.key = key
            self.left_elements = left_elements
            self.right_elements = right_elements
            Bimodule.TensorElement.__init__(self, module, len(left_elements), len(right_elements),
                                            {self: module.ring.one()})

        def __mul__(self, other: AMinus.Generator) -> Bimodule.TensorElement:
            return Bimodule.TensorGenerator(self.module, self.key, self.left_elements, self.right_elements + (other,))

        @multimethod
        def __rmul__(self, other: AMinus.Generator) -> Bimodule.TensorElement:
            return Bimodule.TensorGenerator(self.module, self.key, (other,) + self.left_elements, self.right_elements)

        @multimethod
        def __rmul__(self, other: Z2Polynomial) -> Bimodule.TensorElement:
            return Bimodule.TensorElement.__rmul__(self, other)

        @multimethod
        def __eq__(self, other: Bimodule.TensorGenerator):
            return self.module == other.module and \
                   self.key == other.key and \
                   self.left_elements == other.left_elements and \
                   self.right_elements == other.right_elements

        @multimethod
        def __eq__(self, other):
            return Bimodule.TensorElement.__eq__(self, other)

        def __str__(self):
            return str((self.left_elements, self.key, self.right_elements))

        def __repr__(self):
            return str((self.left_elements, self.key, self.right_elements))

        def __hash__(self):
            return hash((self.left_elements, self.key, self.right_elements))


class TypeDA(Bimodule):
    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 left_scalar_action: Z2PolynomialRing.Map, right_scalar_action: Z2PolynomialRing.Map):
        super().__init__(ring, left_algebra, right_algebra, left_scalar_action, right_scalar_action)

    def to_agraph(self, idempotents=True) -> AGraph:
        graph = AGraph(strict=False, directed=True)
        for generator in self.generators:
            graph.add_node(generator.key,
                           shape='box',
                           fontname='Arial')
        for gen_in, elt_out in self.structure_maps.items():
            for gen_out, c_out in elt_out.coefficients.items():
                graph.add_edge(gen_in.key, gen_out.key,
                               label=str((gen_out.left_elements, c_out, gen_in.right_elements)),
                               dir='forward',
                               color='black',
                               fontname='Arial')
        graph.layout('dot')
        return graph


    # def reduce(self):
    #     reducible_edge = self.reducible_edge()
    #     while reducible_edge is not None:
    #         self.reduce_edge(*reducible_edge)
    #
    # def reducible_edge(self) -> Optional[Tuple]:
    #     pass
    #
    # def reduce_edge(self, source: Bimodule.Generator, target: Bimodule.Generator):
    #     pass

    # def tensor(self, other: TypeDA) -> TypeDA:
    #     assert self.right_algebra.ss == other.left_algebra.ss
    #
    #     generators = [Bimodule.Generator((xm, xn), xm.left_idempotent, xn.right_idempotent)
    #                   for xm in self.graph.nodes for xn in other.graph.nodes
    #                   if xm.right_idempotent == xn.left_idempotent]
    #
    #     in1, in2 = self.ring.tensor_inclusions(other.ring)
    #
    #     maps = set()
    #     for x in generators:
    #         xm = x.key[0]
    #         xn = x.key[1]
    #         for ym in self.graph[xm]:
    #             for i in self.graph[xm][ym]:
    #                 delta_1 = self.graph[xm][ym][i]
    #                 if len(delta_1['left']) > 1:
    #                     continue
    #                 for yn in other.graph[xn]:
    #                     if ym.right_idempotent != yn.left_idempotent:
    #                         continue
    #                     for j in other.graph[xn][yn]:
    #                         delta_n = other.graph[xn][yn][j]
    #                         if delta_1['right'] != delta_n['left']:
    #                             continue
    #                         maps.add(
    #                             Bimodule.Edge(
    #                                 Bimodule.Generator((xm, xn), xm.left_idempotent, xn.right_idempotent),
    #                                 Bimodule.Generator((ym, yn), ym.left_idempotent, yn.right_idempotent),
    #                                 in1.apply(delta_1['c']) * in2.apply(delta_n['c']), delta_1['left'], delta_n['right']))
    #
    #     return TypeDA(in1.target, self.left_algebra, other.right_algebra, generators, maps)
