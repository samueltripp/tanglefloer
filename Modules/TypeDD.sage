class TypeDD:
	# gens - a list of generators
	# edges - a list of Edges
	def __init__(self, gens, edges):
		self.gens = gens
		self.maps = maps
	
	# other is a TypeAA
	def tensor(self, other):
		# TODO
		break
		
	def reduce(self):
		# TODO
		break
		
	# represents the action of some delta_1 on some pair of generators
	class Edge:
		# source, target - elements of gens
		# l_coefficient - element of A
		# r_coefficient - element of B
		# m_cofficient - element of k
		# 
		# condition: delta_1(source) = l_coefficient (X) r_coefficient (X) (m_coefficient * target)
		def __init__(self, source, target, l_coefficient, r_coefficient, m_coefficient)
			self.source = source
			self.target = target
			self.l_coefficient = l_coefficient
			self.r_coefficient = r_coefficient
			self.m_coefficient = m_coefficient
