class TypeDA:
	# gens - a list of generators
	# edges - a list of Edges
	def __init__(self, gens, edges):
		self.gens = gens
		self.maps = maps
	
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
		# r_multipliers - a tuple of j-1 elements of B
		# m_cofficient - element of k
		# l_coefficient - element of A
		# 
		# condition: delta_1^j(source (X) r_multipliers) = l_coefficient (X) (m_coefficient * target)
		def __init__(self, source, target, r_multipliers, m_coefficient, l_coefficient)
			self.source = source
			self.target = target
			self.r_multipliers = r_multipliers
			self.l_coefficient = l_coefficient
			self.m_coefficient = m_coefficient
