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
                 left_scalar_action: Z2PolynomialRing.Map, right_scalar_action: Z2PolynomialRing.Map,
                 generators=None, structure_maps=None):
        self.ring = ring
        self.left_algebra = left_algebra
        self.right_algebra = right_algebra
        self.left_scalar_action = left_scalar_action
        self.right_scalar_action = right_scalar_action
        self.generators = generators or set()  # {Bimodule.Generator}
        self.structure_maps = structure_maps or {}  # {(int, int): {Bimodule.TensorGenerator: Bimodule.TensorElement}}

    def __repr__(self) -> str:
        return str(self.__dict__)

    def add_generator(self, generator: Bimodule.Generator) -> None:
        self.generators.add(generator)

    def add_structure_map(self, input: Bimodule.Generator, output: Bimodule.Element) -> None:
        pass

    # returns the zero element of A^(x)i (x) M (x) A^(x)j
    def zero(self, i=0, j=0) -> Bimodule.Element:
        return Bimodule.Element(self, i, j, {})

    # represents an element of A^(x)i (x) M (x) A^(x)j
    class Element:
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

        def __add__(self, other: Bimodule.Element) -> Bimodule.Element:
            assert self.i == other.i and self.j == other.j
            new_coefficients = dict(self.coefficients)
            for g in other.coefficients:
                if g in self.coefficients:
                    new_coefficients[g] = self.coefficients[g] + other.coefficients[g]
                else:
                    new_coefficients[g] = other.coefficients[g]
            return Bimodule.Element(self.module, self.i, self.j, new_coefficients)

        @multimethod
        def __rmul__(self, other: Z2Polynomial) -> Bimodule.Element:
            new_coefficients = dict(self.coefficients)
            for g in new_coefficients:
                new_coefficients[g] = other * new_coefficients[g]
            return Bimodule.Element(self.module, self.i, self.j, new_coefficients)

        @multimethod
        def __rmul__(self, other: AMinus.Element) -> Bimodule.Element:
            out = self.module.zero(self.i + 1, self.j)

            for g1, c1 in other.coefficients.items():
                for g2, c2 in self.coefficients.items():
                    out += (self.module.left_scalar_action.apply(c1) * c2) * (g1 * g2).to_element()

            return out

        def __mul__(self, other: AMinus.Element) -> Bimodule.Element:
            out = self.module.zero(self.i, self.j + 1)

            for g1, c1 in self.coefficients.items():
                for g2, c2 in other.coefficients.items():
                    out += (c1 * self.module.right_scalar_action.apply(c2)) * (g1 * g2)

            return out

        def __eq__(self, other: Bimodule.Element) -> bool:
            return self.module == other.module and self.coefficients == other.coefficients

        def __hash__(self):
            return hash((self.module, self.coefficients))

        def __repr__(self) -> str:
            return str(self.coefficients)

    # represents a generator of A^(x)i (x) M (x) A^(x)j
    class Generator:
        def __init__(self, module: Bimodule, key, left_idempotent: AMinus.Generator, right_idempotent: AMinus.Generator,
                     left: Tuple[AMinus.Generator, ...] = None, right: Tuple[AMinus.Generator, ...] = None):
            self.module = module
            self.key = key
            self.left_idempotent = left_idempotent  # the left idempotent of m
            self.right_idempotent = right_idempotent  # the right idempotent of m
            self.left = left or tuple()
            self.right = right or tuple()

        def to_element(self) -> Bimodule.Element:
            return Bimodule.Element(self.module, len(self.left), len(self.right), {self: self.module.ring.one()})

        def __mul__(self, other: AMinus.Generator) -> Bimodule.Generator:
            return Bimodule.Generator(self.module, self.key, self.left_idempotent, self.right_idempotent,
                                      self.left, self.right + (other,))

        def __rmul__(self, other: AMinus.Generator) -> Bimodule.Generator:
            return Bimodule.Generator(self.module, self.key, self.left_idempotent, self.right_idempotent,
                                      (other,) + self.left, self.right)

        def __str__(self):
            return str((self.left, self.key, self.right))

        def __repr__(self):
            return str((self.left, self.key, self.right))

        def __eq__(self, other: Bimodule.Generator):
            return self.module == other.module and \
                   self.key == other.key and \
                   self.left == other.left and \
                   self.right == other.right

        def __hash__(self):
            return hash((self.module, self.key, self.left, self.right))


class TypeDA(Bimodule):
    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 left_scalar_action: Z2PolynomialRing.Map, right_scalar_action: Z2PolynomialRing.Map,
                 generators=None, structure_maps=None):
        super().__init__(ring, left_algebra, right_algebra, left_scalar_action, right_scalar_action,
                         generators=None, structure_maps=None)

    def add_structure_map(self, input: Bimodule.Generator, output: Bimodule.Element) -> None:
        assert len(input.left) == output.j == 0
        size = (output.i, len(input.right))
        if size not in self.structure_maps:
            self.structure_maps[size] = {}
        if input not in self.structure_maps[size]:
            self.structure_maps[size][input] = self.zero(output.i, 0)
        self.structure_maps[size][input] += output

    def to_agraph(self, idempotents=True) -> AGraph:
        graph = AGraph(strict=False, directed=True)
        for generator in self.generators:
            graph.add_node(generator.key,
                           shape='box',
                           fontname='Arial')
        for size in self.structure_maps.keys():
            for gen_in, elt_out in self.structure_maps[size].items():
                if not idempotents and len(gen_in.right) == 1 and gen_in.right[0].is_idempotent():
                    continue
                for gen_out, c_out in elt_out.coefficients.items():
                    graph.add_edge(gen_in.key, gen_out.key,
                                   label=str((gen_out.left, c_out, gen_in.right)),
                                   dir='forward',
                                   color='blue' if len(gen_in.right) > 0 else 'black',
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

    def tensor(self, other: TypeDA) -> TypeDA:
        assert self.right_algebra.ss == other.left_algebra.ss

        in_m, in_n = self.ring.tensor_inclusions(other.ring)

        out = TypeDA(in_m.target, self.left_algebra, other.right_algebra,
                     in_m.compose(self.left_scalar_action), in_n.compose(other.right_scalar_action))

        for xm in self.generators:
            for xn in other.generators:
                if xm.right_idempotent == xn.left_idempotent:
                    out.add_generator(Bimodule.Generator(out, (xm.key, xn.key), xm.left_idempotent, xn.right_idempotent))

        for size_m in self.structure_maps.keys():
            for xm, yms in self.structure_maps[size_m].items():
                for ym, cm in yms.coefficients.items():
                    if len(ym.left) != 1:
                        continue
                    for size_n in other.structure_maps.keys():
                        for xn, yns in other.structure_maps[size_n].items():
                            if xm.right_idempotent != xn.left_idempotent:
                                continue
                            for yn, cn in yns.coefficients.items():
                                if ym.right_idempotent != yn.left_idempotent:
                                    continue
                                if xm.right != yn.left:
                                    continue
                                x = Bimodule.Generator(out, (xm.key, xn.key), xm.left_idempotent, xn.right_idempotent,
                                                       xm.left, xn.right)
                                y = in_m.apply(cm) * in_n.apply(cn) * \
                                    Bimodule.Generator(out, (ym.key, yn.key), ym.left_idempotent, yn.right_idempotent,
                                                       ym.left, yn.right).to_element()
                                out.add_structure_map(x, y)

        return out
