class TypeDD:
    # gens - a list of generators
    # edges_out - a dictionary keyed by source generators
    # edges_in - a dictionary keyed by target generators
    # deleted edge data...it is stored in the gens
    # dictionary values are lists of edges (each edge is a 5-tuple (source, target, a, b, m))
    def __init__(self, gens, edges_out, edges_in):
        self.gens = gens
        self.edges_out = edges_out
        self.edges_in = edges_in

    # other is a TypeAA
    def tensor(self, other):
        # TODO
        break

    # test if edge is empty
    # maybe this should go somewhere else ?
    # input is an edge
    def is_empty_edge(e):
        assertion_1 = (e.source.idempotent_left == e.a_coefficient)
        assertion_2 = (e.source.idempotent_right == e.a_coefficient)
        assertion_3 = (e.m_coefficient == 1) # 1 in kk
        return assertion_1 and assertion_2 and assertion_3

    # input: TypeDD
    # len(output) in {0,1}
    def get_reducible_edge(self):
        for x in self.gens: # x is a generator
            for e in x.edges_out:
                if is_empty_edge(e):
                    # now check to see if there are other edges between x and target of e
                    y = e[1] # target of e
                    # check size of list of edges from x to y is 1
                    assertion_1 = len(get_edges_between(self, x, y)) == 1
                    # check that there are no edges from y to x
                    assertion_2 = len(get_edges_between(self, y, x)) == 0
                    if assertion_1 and assertion_2:
                        return [e]
        return []

    # get edges z->w
    # self: TypeDD
    # z: source generator
    # w: target generator
    # returns list of edges z->w
    def get_edges_between(self, z, w):
        edges_out_z = self.edges_out(z)
        edges_in_w = self.edges_in(w)
        edges = []
        for edge in edges_out_z:
            if edge in edges_in_w:
                edges.append(edge)
        return edges

    # input is a TypeDD
    # copy generators and update
    # or we could build a new TypeDD completely...
    def reduce(self):
        edge_list = get_reducible_edge(self)
        while len(edge_list) != 0:
            edge = edge_list[0] # take the first reducible edge and reduce it
            x = edge.source
            y = edge.target
            # new list of generators without reference to self.gens
            new_gens = []
            for v in self.gens:
                if v != x and v != y:
                    new_gens.append(v)
            # new edges_out edges_in dictionaries each with an empty list of edges
            new_edges_out = {}
            new_edges_in = {}
            for v in new_gens:
                new_edges_out[v] = []
                new_edges_in[v] = []
            for z in new_gens:
                # compute all targets of z
                z_targets = [ edge[1] for edge in z.edges_out ]
                for w in new_gens:
                    # get edge(s) z->w
                    edges_z_to_w = get_edges_between(self, z, w)
                    z_points_to_y = y in z_targets
                    # compute all sources of w
                    w_sources = [ edge[0] for edge in w.edges_in ]
                    x_points_to_w = x in w_sources
                    # if z points to y and x points to x then update edge from z to w
                    if z_points_to_y and x_points_to_w: # otherwise don't update anything
                        # TODO: a, b, m determined by MATH
                        new_edge = Edge(z, w, a, b, m)
                        new_edges_out[z].append(new_edge)
                        new_edges_in[w].append(new_edge)
                    else: # just copy over from previous dictionaries
                        new_edges_out[z] += edges_z_to_w
                        new_edges_in[w] += edges_z_to_w

            # delete source and target generators of edge to be reduced

            # update self
            self.gens = new_gens
            self.edges_out = new_edges_out
            self.edges_in = new_edges_in
            # update (reducible) edge_list
            edge_list = get_reducible_edge(self)
        assert len(edge_list) == 0 # assert that there are no more edges to reduce
        return self

    class Generator:
        # in/out edges
        # left/right idempotents
        def __init__(self, edges_in, edges_out, idempotent_left, idempotent_right):
            self.edges_in = edges_in
            self.edges_out = edges_out
            self.idempotent_left = idempotent_left
            self.idempotent_right = idempotent_right

    # represents the action of some delta_1 on some pair of generators
    class Edge:
        # source, target - elements of gens
        # a_coefficient - element of A
        # m_cofficient - element of k
        # b_coefficient - element of B
        # 
        # condition: delta_1(source) = a_coefficient (X) (m_coefficient * target) (X) b_coefficient
        def __init__(self, source, target, a_coefficient, b_coefficient, m_coefficient):
            self.source = source
            self.target = target
            self.a_coefficient = a_coefficient
            self.b_coefficient = b_coefficient
            self.m_coefficient = m_coefficient
