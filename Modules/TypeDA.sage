class TypeDA:
	# gens - a list of generators
	# edges - a list of Edges
	def __init__(self, gens, edges):
		self.gens = gens
		self.edges = edges
	
	# two cases: other is a TypeDA, or other is a TypeDD
	def tensor(self, other):
		# TODO
		break
		
	def reduce(self):
		# TODO
		break
		
	# represents the action of some delta_1^j on some pair of generators
	class Edge:
		# source, target - elements of gens
		# b_multipliers - a tuple of j-1 elements of B
		# m_cofficient - element of k
		# a_coefficient - element of A
		# 
		# condition: delta_1^j(source (X) b_multipliers) = a_coefficient (X) (m_coefficient * target)
		def __init__(self, source, target, b_multipliers, a_coefficient, m_coefficient)
			self.source = source
			self.target = target
			self.b_multipliers = b_multipliers
			self.a_coefficient = a_coefficient
			self.m_coefficient = m_coefficient
