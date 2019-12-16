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

    def add_structure_map(self, input: Bimodule.TensorGenerator, output: Bimodule.TensorElement) -> None:
        pass

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
            if self.i != other.i or self.j != other.j:
                print('hi')
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
            return self.module == other.module and self.coefficients == other.coefficients

        def __hash__(self):
            return hash((self.module, self.coefficients))

        def __repr__(self) -> str:
            return str(self.coefficients)

    # represents a generator of A^(x)i (x) M (x) A^(x)j
    class TensorGenerator(TensorElement):
        def __init__(self, module: Bimodule, generator: Bimodule.Generator,
                     left_generators: Tuple[AMinus.Generator, ...], right_generators: Tuple[AMinus.Generator, ...]):
            self.module = module
            self.generator = generator
            self.left_generators = left_generators
            self.right_generators = right_generators
            Bimodule.TensorElement.__init__(self, module, len(left_generators), len(right_generators),
                                            {self: module.ring.one()})

        def __mul__(self, other: AMinus.Generator) -> Bimodule.TensorElement:
            return Bimodule.TensorGenerator(self.module, self.generator, self.left_generators,
                                            self.right_generators + (other,))

        @multimethod
        def __rmul__(self, other: AMinus.Generator) -> Bimodule.TensorElement:
            return Bimodule.TensorGenerator(self.module, self.generator, (other,) + self.left_generators,
                                            self.right_generators)

        @multimethod
        def __rmul__(self, other: Z2Polynomial) -> Bimodule.TensorElement:
            return Bimodule.TensorElement.__rmul__(self, other)

        @multimethod
        def __eq__(self, other: Bimodule.TensorGenerator):
            return self.module == other.module and \
                   self.generator == other.generator and \
                   self.left_generators == other.left_generators and \
                   self.right_generators == other.right_generators

        @multimethod
        def __eq__(self, other):
            return Bimodule.TensorElement.__eq__(self, other)

        def __str__(self):
            return str((self.left_generators, self.generator, self.right_generators))

        def __repr__(self):
            return str((self.left_generators, self.generator, self.right_generators))

        def __hash__(self):
            return hash((self.module, self.generator, self.left_generators, self.right_generators))

    class Element(TensorElement):
        def __init__(self, module: Bimodule, i, j, coefficients=None):
            super().__init__(module, i, j, coefficients)

    class Generator(TensorGenerator, Element):
        def __init__(self, module: Bimodule, key,
                     left_idempotent: AMinus.Generator, right_idempotent: AMinus.Generator):
            self.module = module
            self.key = key
            self.left_idempotent = left_idempotent
            self.right_idempotent = right_idempotent
            Bimodule.TensorGenerator.__init__(self, module, self, tuple(), tuple())
            Bimodule.Element.__init__(self, module, 0, 0, {self: module.ring.one()})

        def __repr__(self):
            return str(self.key)

        def __str__(self):
            return str(self.key)

        @multimethod
        def __eq__(self, other) -> bool:
            return Bimodule.TensorElement.__eq__(self, other)

        @multimethod
        def __eq__(self, other: Bimodule.Generator) -> bool:
            return self.module == other.module and \
                   self.key == other.key and \
                   self.left_idempotent == other.left_idempotent and \
                   self.right_idempotent == other.right_idempotent

        def __hash__(self):
            return hash((self.module, self.key, self.left_idempotent, self.right_idempotent))


class TypeDA(Bimodule):
    def __init__(self, ring: Z2PolynomialRing, left_algebra: AMinus, right_algebra: AMinus,
                 left_scalar_action: Z2PolynomialRing.Map, right_scalar_action: Z2PolynomialRing.Map,
                 generators=None, structure_maps=None):
        super().__init__(ring, left_algebra, right_algebra, left_scalar_action, right_scalar_action,
                         generators=None, structure_maps=None)

    def add_structure_map(self, input: Bimodule.TensorGenerator, output: Bimodule.TensorElement) -> None:
        assert input.i == output.j == 0
        size = (output.i, input.j)
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
            for tgen_in, telt_out in self.structure_maps[size].items():
                if not idempotents and len(tgen_in.right_generators) == 1 and tgen_in.right_generators[0].is_idempotent():
                    continue
                for tgen_out, c_out in telt_out.coefficients.items():
                    graph.add_edge(tgen_in.generator.key, tgen_out.generator.key,
                                   label=str((tgen_out.left_generators, c_out, tgen_in.right_generators)),
                                   dir='forward',
                                   color='blue' if len(tgen_in.right_generators) > 0 else 'black',
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
                    out.add_generator(Bimodule.Generator(out, (xm, xn), xm.left_idempotent, xn.right_idempotent))

        for size_m in self.structure_maps.keys():
            for txm, tyms in self.structure_maps[size_m].items():
                xm = txm.generator
                for tym, cm in tyms.coefficients.items():
                    ym = tym.generator
                    if len(tym.left_generators) != 1:
                        continue
                    for size_n in other.structure_maps.keys():
                        for txn, tyns in other.structure_maps[size_n].items():
                            xn = txn.generator
                            if xm.right_idempotent != xn.left_idempotent:
                                continue
                            for tyn, cn in tyns.coefficients.items():
                                yn = tyn.generator
                                if ym.right_idempotent != yn.left_idempotent:
                                    continue
                                if xm.right_generators != yn.left_generators:
                                    continue
                                x = Bimodule.Generator(out, (xm, xn), xm.left_idempotent, xn.right_idempotent)
                                y = Bimodule.Generator(out, (ym, yn), ym.left_idempotent, yn.right_idempotent)
                                out.add_structure_map(
                                    Bimodule.TensorGenerator(out, x, txm.left_generators, txn.right_generators),
                                    in_m.apply(cm) * in_n.apply(cn) *
                                    Bimodule.TensorGenerator(out, y, tym.left_generators, tyn.right_generators))

        return out
