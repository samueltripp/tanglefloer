from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Optional, List

import networkx as nx
from frozendict import frozendict
from multimethod import multimethod
from networkx import DiGraph, MultiDiGraph
from pygraphviz import AGraph

from Functions.Functions import simplify_coefficients
from SignAlgebra.AMinus import AMinus
from SignAlgebra.TensorAlgebra import TensorAlgebra
from SignAlgebra.Z2PolynomialRing import Z2PolynomialRing, Z2Polynomial


# Base class for chain complexes and Type D, A, DD, AA, DA, and AD structures
class Module(ABC):
    def __init__(self, ring: Z2PolynomialRing, left_algebra: Optional[AMinus], right_algebra: Optional[AMinus],
                 left_scalar_action: Optional[Z2PolynomialRing.Map],
                 right_scalar_action: Optional[Z2PolynomialRing.Map],
                 graph: MultiDiGraph = None):
        self.ring = ring
        self.left_algebra = left_algebra
        self.right_algebra = right_algebra
        self.left_scalar_action = left_scalar_action
        self.right_scalar_action = right_scalar_action
        self.graph = graph or MultiDiGraph()
        if left_scalar_action is not None:
            self.left_tensor_algebra = TensorAlgebra(left_algebra)
            assert left_algebra.ring == left_scalar_action.source, "left scalar action has wrong source"
            assert ring == left_scalar_action.target, "left scalar action has wrong target"
        if right_scalar_action is not None:
            self.right_tensor_algebra = TensorAlgebra(right_algebra)
            assert right_algebra.ring == right_scalar_action.source, "right scalar action has wrong source"
            assert ring == right_scalar_action.target, "right scalar action has wrong target"

    # add the structure map (input |-> output) to this module
    def add_structure_map(self, input: Module.TensorGenerator, output: Module.TensorElement) -> None:
        assert self.valid_input_gen(input)
        x = input.get_module_generator()
        for gen_out, c_out in output.coefficients.items():
            assert self.valid_output_gen(gen_out)
            y = gen_out.get_module_generator()
            self.add_edge(x, y, (gen_out.left, input.right), c_out)

    @staticmethod
    def valid_input_gen(g):
        pass

    @staticmethod
    def valid_output_gen(g):
        pass

    def __repr__(self) -> str:
        return str(self.__dict__)

    # add the given generator to this module
    def add_generator(self, generator: Module.TensorGenerator, grading: (int, int)) -> None:
        self.graph.add_node(generator, grading=grading)

    def add_edge(self, x, y, k, c):
        current = self.graph.get_edge_data(x, y, key=k)
        if current is None:
            self.graph.add_edge(x, y, key=k, c=self.ring.zero())
            current = self.graph.get_edge_data(x, y, key=k)
        current['c'] += c
        if current['c'] == self.ring.zero():
            self.graph.remove_edge(x, y, key=k)

    def reduce(self) -> Module:
        components = self.decomposed()
        if len(components) == 0:
            return self
        for component in components:
            component.reduce_component()
        return self.direct_sum(components)

    def reduce_component(self) -> Module:
        reducible_edge = self.get_reducible_edge()
        while reducible_edge is not None:
            self.reduce_edge(*reducible_edge)
            reducible_edge = self.get_reducible_edge()
        return self

    def get_reducible_edge(self):
        for x in self.graph:
            for y in self.graph[x]:
                if len(self.graph[x][y]) == 1:
                    k, d = list(self.graph[x][y].items())[0]
                    left = k[0]
                    right = k[1]
                    if self.edge_is_reducible(left, d['c'], right):
                        return x, y, k, d

    @staticmethod
    @abstractmethod
    def edge_is_reducible(left, c, right):
        pass

    @abstractmethod
    def reduce_edge(self, x, y, k, d):
        pass

    def identify_variables(self, var1, var2):
        assert var1 in self.ring.variables and var2 in self.ring.variables
        r_merged = Z2PolynomialRing([v for v in self.ring.variables if v != var2])
        f_merge = Z2PolynomialRing.Map(self.ring, r_merged,
                                       {v: (v if v != var2 else var1) for v in self.ring.variables})

        left_scalar_action_merged = f_merge.compose(self.left_scalar_action) if self.left_scalar_action else None
        right_scalar_action_merged = f_merge.compose(self.right_scalar_action) if self.right_scalar_action else None

        subclass = type(self)
        out = subclass(r_merged, self.left_algebra, self.right_algebra,
                       left_scalar_action_merged, right_scalar_action_merged)

        for x, x_data in self.graph.nodes(data=True):
            out.add_generator(Module.change_module_of_generator(x, out), grading=x_data['grading'])
        for x in self.graph.nodes():
            for _, y, k, xyk_data in self.graph.out_edges(x, keys=True, data=True):
                c = xyk_data['c']
                c = f_merge.apply(c)
                out.add_edge(Module.change_module_of_generator(x, out), Module.change_module_of_generator(y, out), k, c)
        return out

    def copy(self):
        subclass = type(self)
        return subclass(self.ring,
                        self.left_algebra, self.right_algebra,
                        self.left_scalar_action, self.right_scalar_action,
                        self.graph.copy())

    # mostly a helper method for identify_variables()
    @staticmethod
    def change_module_of_generator(x, m):
        return Module.TensorGenerator(m, x.key, x.left_idempotent, x.right_idempotent, x.left, x.right)

    # shift gradings
    def __getitem__(self, item):
        maslov_shift, ta_shift = item
        m = self.copy()
        for x, x_data in self.graph.nodes(data=True):
            maslov, ta = x_data['grading']
            m.graph.add_node(x, grading=(maslov + maslov_shift, ta + ta_shift))
        return m

    def is_isomorphic_to(self, other: Module):
        if (self.ring != other.ring) \
                or (self.left_algebra != other.left_algebra) \
                or (self.right_algebra != other.right_algebra) \
                or (self.left_scalar_action != other.left_scalar_action) \
                or (self.right_scalar_action != other.right_scalar_action):
            return False
        node_match = lambda x_data, y_data: x_data == y_data
        edge_match = lambda e1_data, e2_data: e1_data == e2_data
        return nx.is_isomorphic(self.graph, other.graph, node_match=node_match, edge_match=edge_match)

    # if self is homotopic to C (+) C[1,2] for some C, return C, else return None
    def halve(self):
        components = self.decomposed()
        component_matches = DiGraph()
        for i in range(len(components)):
            component_matches.add_node(i)

        for i, comp1 in enumerate(components):
            for j, comp2 in enumerate(components):
                if comp1.is_isomorphic_to(comp2[1, 2]):
                    component_matches.add_edge(i, j)

        matching = nx.max_weight_matching(component_matches.to_undirected(), maxcardinality=True)
        if len(matching) < len(components) / 2:
            return None
        matching_sources = [i if component_matches.has_edge(i, j) else j for (i, j) in matching]
        return Module.direct_sum([components[i] for i in matching_sources])

    def decomposed(self):
        subclass = type(self)
        return [subclass(self.ring, self.left_algebra, self.right_algebra,
                         self.left_scalar_action, self.right_scalar_action,
                         MultiDiGraph(self.graph.subgraph(component)))
                for component in nx.weakly_connected_components(self.graph)]

    @staticmethod
    def direct_sum(modules: List):
        subclass = type(modules[0])
        new_graph = nx.union_all([m.graph for m in modules])
        return subclass(modules[0].ring, modules[0].left_algebra, modules[0].right_algebra,
                        modules[0].left_scalar_action, modules[0].right_scalar_action, new_graph)

    def zero(self) -> Module.TensorElement:
        return Module.TensorElement(self)

    # turns this bimodule into a graphviz-compatible format
    def to_agraph(self, idempotents=True) -> AGraph:
        graph = AGraph(strict=False, directed=True)
        for generator, data in self.graph.nodes(data=True):
            graph.add_node(str(generator.key) + '[' + str(data['grading'])[1:-1] + ']',
                           shape='box',
                           fontname='Arial')
        for x, y, (left, right), d in self.graph.edges(keys=True, data=True):
            c = d['c']
            x_grading = self.graph.nodes(data=True)[x]['grading']
            y_grading = self.graph.nodes(data=True)[y]['grading']
            if not idempotents and self.is_idempotent_edge_data(left, c, right):
                continue
            graph.add_edge(str(x.key) + '[' + str(x_grading)[1:-1] + ']',
                           str(y.key) + '[' + str(y_grading)[1:-1] + ']',
                           label=' ' + str((left, c, right)) + ' ',
                           dir='forward',
                           color=self.edge_color(left, c, right),
                           fontname='Arial',
                           decorate=True)
        graph.layout('dot')
        return graph

    @staticmethod
    @abstractmethod
    def is_idempotent_edge_data(left, c, right):
        pass

    @staticmethod
    @abstractmethod
    def edge_color(left, c, right):
        pass

    class TensorElement:
        # coefficients - {TensorGenerator: Z2Polynomial}
        def __init__(self, module: Module, coefficients=None):
            self.module = module
            if coefficients is None:
                coefficients = {}
            self.coefficients = frozendict(simplify_coefficients(coefficients))

        def __add__(self, other: Module.TensorElement) -> Module.TensorElement:
            new_coefficients = dict(self.coefficients)
            for g in other.coefficients:
                if g in self.coefficients:
                    new_coefficients[g] = self.coefficients[g] + other.coefficients[g]
                else:
                    new_coefficients[g] = other.coefficients[g]
            return Module.TensorElement(self.module, new_coefficients)

        def __rmul__(self, other: Z2Polynomial) -> Module.TensorElement:
            if (self.module.left_algebra is not None) and (other.ring == self.module.left_algebra.ring) \
                    and (self.module.left_scalar_action is not None):
                other = self.module.left_scalar_action.apply(other)
            elif (self.module.right_algebra is not None) and (other.ring == self.module.right_algebra.ring) \
                    and (self.module.right_scalar_action is not None):
                other = self.module.right_scalar_action.apply(other)
            if other.ring == self.module.ring:
                new_coefficients = {}
                for g, c in self.coefficients.items():
                    new_coefficients[g] = other * c
                return Module.TensorElement(self.module, new_coefficients)
            else:
                raise Exception("no action by the given polynomial")

        @multimethod
        def __rpow__(self, other: AMinus.Generator) -> Module.TensorElement:
            return other.to_element() ** self

        # the tensor product A (x) (A^(x)i (x) M (x) A^(x)j) -> A^(x)i+1 (x) M (x) A^(x)j
        @multimethod
        def __rpow__(self, other: AMinus.Element) -> Module.TensorElement:
            out = self.module.zero()

            for g1, c1 in other.coefficients.items():
                for g2, c2 in self.coefficients.items():
                    if g1.right_idempotent() == g2.leftmost_idempotent():
                        if self.module.left_scalar_action is not None:
                            out += (self.module.left_scalar_action.apply(c1) * c2) * (g1 ** g2)
                        else:
                            raise Exception('no left scalar action')

            return out

        @multimethod
        def __rpow__(self, other: TensorAlgebra.Generator):
            return other.to_element() ** self

        @multimethod
        def __rpow__(self, other: TensorAlgebra.Element):
            out = self.module.zero()

            for g1, c1 in other.coefficients.items():
                for g2, c2 in self.coefficients.items():
                    if g1.right_idempotent() == g2.leftmost_idempotent():
                        if self.module.left_scalar_action is not None:
                            out += (self.module.left_scalar_action.apply(c1) * c2) * (g1 ** g2)
                        else:
                            raise Exception('no left scalar action')

            return out

        # the tensor product (A^(x)i (x) M (x) A^(x)j) (x) A -> A^(x)i (x) M (x) A^(x)j+1
        def __pow__(self, other: AMinus.Element) -> Module.TensorElement:
            out = self.module.zero()

            for g1, c1 in self.coefficients.items():
                for g2, c2 in other.coefficients.items():
                    if g1.rightmost_idempotent() == g2.left_idempotent():
                        if self.module.right_scalar_action is not None:
                            out += (c1 * self.module.right_scalar_action.apply(c2)) * (g1 ** g2)
                        else:
                            raise Exception('no right scalar action')

            return out

        @multimethod
        def __eq__(self, other: Module.TensorElement) -> bool:
            return self.module == other.module and self.coefficients == other.coefficients

        @multimethod
        def __eq__(self, other: Module.TensorGenerator) -> bool:
            return self == other.to_element()

        def __hash__(self):
            return hash((self.module, self.coefficients))

        def __repr__(self) -> str:
            return str(dict(self.coefficients))

    # represents a generator of A^(x)i (x) M (x) A^(x)j as a module over the base ring of M
    class TensorGenerator:
        def __init__(self, module: Module, key,
                     left_idempotent: AMinus.Generator, right_idempotent: AMinus.Generator,
                     left: TensorAlgebra.Generator = None,
                     right: TensorAlgebra.Generator = None):
            self.module = module
            self.key = key
            self.left_idempotent = left_idempotent
            self.right_idempotent = right_idempotent
            self.left = left or module.left_tensor_algebra.one_generator()
            self.right = right or module.right_tensor_algebra.one_generator()

        # converts this generator to an actual element
        def to_element(self) -> Module.TensorElement:
            return Module.TensorElement(self.module, {self: self.module.ring.one()})

        def get_module_generator(self) -> Module.TensorGenerator:
            return Module.TensorGenerator(self.module, self.key, self.left_idempotent, self.right_idempotent)

        def num_left_factors(self):
            return self.left.num_factors()

        def num_right_factors(self):
            return self.right.num_factors()

        def leftmost_idempotent(self) -> AMinus.Generator:
            if self.left.num_factors() == 0:
                return self.left_idempotent
            else:
                return self.left[0].left_idempotent()

        def rightmost_idempotent(self) -> AMinus.Generator:
            if self.right.num_factors() == 0:
                return self.right_idempotent
            else:
                return self.right[-1].right_idempotent()

        # addition
        def __add__(self, other) -> Module.TensorElement:
            return self.to_element() + other

        def __rmul__(self, other: Z2Polynomial) -> Module.TensorElement:
            return other * self.to_element()

        # tensor product
        @multimethod
        def __pow__(self, other: AMinus.Generator) -> Module.TensorGenerator:
            assert self.rightmost_idempotent() == other.left_idempotent()
            return Module.TensorGenerator(self.module, self.key, self.left_idempotent, self.right_idempotent,
                                          self.left, self.right ** other)

        # tensor product
        @multimethod
        def __pow__(self, other: AMinus.Element) -> Module.TensorElement:
            return self.to_element() ** other

        # tensor product
        # assumes other is a tuple of generators or elements
        @multimethod
        def __pow__(self, other: TensorAlgebra.Generator) -> Module.TensorGenerator:
            return Module.TensorGenerator(self.module, self.key, self.left_idempotent, self.right_idempotent, self.left,
                                          self.right ** other)

        @multimethod
        def __pow__(self, other: TensorAlgebra.Element) -> Module.TensorElement:
            out = self.module.zero()
            for gen, coefficient in other.coefficients.items():
                out += coefficient * (self ** gen)
            return out

        # tensor product
        @multimethod
        def __rpow__(self, other: AMinus.Generator) -> Module.TensorGenerator:
            assert other.right_idempotent() == self.leftmost_idempotent()
            return Module.TensorGenerator(self.module, self.key,
                                          self.left_idempotent, self.right_idempotent, other ** self.left, self.right)

        # tensor product
        @multimethod
        def __rpow__(self, other: AMinus.Element) -> Module.TensorElement:
            return other ** self.to_element()

        # tensor product
        # assumes other is a tuple of generators or elements
        @multimethod
        def __rpow__(self, other: TensorAlgebra.Generator) -> Module.TensorGenerator:
            return Module.TensorGenerator(self.module, self.key, self.left_idempotent, self.right_idempotent,
                                          other ** self.left, self.right)

        def __str__(self):
            return str((self.left, self.key, self.right))

        def __repr__(self):
            return str((self.left, self.key, self.right))

        @multimethod
        def __eq__(self, other: Module.TensorGenerator):
            return self.module == other.module and \
                   self.key == other.key and \
                   self.left_idempotent == other.left_idempotent and \
                   self.right_idempotent == other.right_idempotent and \
                   self.left == other.left and \
                   self.right == other.right

        @multimethod
        def __eq__(self, other):
            return self.to_element() == other

        def __hash__(self):
            return hash((self.left, self.module, self.key,
                         self.left_idempotent, self.right_idempotent, self.right))
