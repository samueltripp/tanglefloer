import numpy as np

class TypeDD:
	# gens - a list of generators
    # edges_out - a dictionary keyed by source generators
    # edges_in - a dictionary keyed by target generators
    # dictionary values are lists of edges (each edge is a 5-tuple (source, target, a, b, m))
	def __init__(self, gens, edges_out, edges_in):
		self.gens = gens
		
		#self.edges = dictEdges(edges)
		self.edges_out = edges_out
		self.edges_in = edges_in
	
	# other is a TypeAA (only option for TypeDD)
	def tensor(self, other):
		assert self.rightalgebra = other.leftalgebra, #pseudo - where are the left and right algebras stored?
			"Error: Algebras not compatible"
		MNgens = np.zeros(length(self.gens),length(other.gens)) # I think a list will be more efficient than a matrix, since it will probably be quite sparse, but we should discuss
		for i in length(self.gens):
			for j in length(other.gens):
				if self.gens[i][2] == other.gens[j][1]: #this assumes generators are a tuple (x, eL, eR)
					MNgens[i][j] = [self.gens[i][0]+other.gens[j][0],self.gens[i][1],other.gens[j][2]] #(x tensor y, xeL, yeR)
					# above "+" assumes the x generator is stored as a list, which I think it no longer is
					#(but also I think the + operator might have been overwritten for generators, so might work anyway. Check this.)

		# TODO: Change the below to match the matrix MNgens above (used to be a list)
		for xy in MNgens: # look at each (x tensor y) generator
			for e in other.edges[xy[4]]: # looking at each edge coming out of y (stored as a list in dict with key of source vertex)
				# starting to like matrices much better... should go back and change that #TODO
				break

	def reduce(self):
		# TODO
		break

	#creates a dictionary of the edges, keyed by source vertex. The values are a list of edges.
	def dictEdges(self,edges): # matrix instead??
		edgeDict = dict()
		for e in edges:
			if e.source in edgeDict:
				edgeDict[e.source].append(e)
				# I'm not convinced the above works as I intend - need to test small example
			else:
				edgeDict[e.source] = [e]
		return edgeDict

		
    # test if edge is empty
    # maybe this should go somewhere else ?
    def is_empty_edge(self):
        assertion_1 = (self.source.idempotent_left == self.a_coefficient)
        assertion_2 = (self.source.idempotent_right == self.a_coefficient)
        assertion_3 = (self.m_coefficient == 1) # 1 in kk
        if assertion_1 and assertion_2 and assertion_3:
            return True
        else:
            return False

	def get_reducible_edge(self):
        for x in self.gens:
            for e in x.edges_out:
                if is_empty_edge(e):
                    return [e]
        return []

        # for i in range(0, len(gens)):
        #     edge = edges

        # edges_reduced = {} # initialize empty dictionary of new edges
        # for i in range(0, len(edges)): # loop over edges
        #     if is_empty_edge(edges[i])
        #     # for each edge loop over pairs of generators
        #     for j in range(0, len(gens)):
        #         for k in range(0, len(gens)):
        #             # TODO: at least one line of math goes here :)
        # self.edges = edges_reduced
		# break

    class Generator:
    # in/out edges
    # left/right idempotents
    def __init__(self, edges_in, edges_out, idempotent_left, idempotent_right):
        self.edges_in = edges_in
        self.edges_out = edges_out
        self.idempotent_left = idempotent_left
        self.idempotent_right = idempotent_right

>>>>>>> upstream/master
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