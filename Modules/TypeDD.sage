class TypeDD:
	# gens - a list of generators
	# edges - a list of Edges
	def __init__(self, gens, edges):
		self.gens = gens
		self.edges = edges
	
	# other is a TypeAA
	def tensor(self, other):
		# TODO
		break
		
	def reduce(self):
        gens = self.gens # get generators
        edges = self.edges # get edges
        edges_reduced = {} # initialize empty dictionary of new edges
        for i in range(0, len(edges)): # loop over edges
            edge = edges[i]
            # for each edge loop over pairs of generators
            for j in range(0, len(gens)):
                for k in range(0, len(gens)):
                    # TODO: at least one line of math goes here :)
        self.edges = edges_reduced
		break
		
	# represents the action of some delta_1 on some pair of generators
	class Edge:
		# source, target - elements of gens
		# a_coefficient - element of A
		# m_cofficient - element of k
		# b_coefficient - element of B
		# 
		# condition: delta_1(source) = a_coefficient (X) (m_coefficient * target) (X) b_coefficient
		def __init__(self, source, target, a_coefficient, b_coefficient, m_coefficient)
			self.source = source
			self.target = target
			self.a_coefficient = a_coefficient
			self.b_coefficient = b_coefficient
			self.m_coefficient = m_coefficient
