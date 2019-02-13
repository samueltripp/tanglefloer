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
		# a_coefficient - element of A
		# b_coefficient - element of B
		# m_cofficient - element of k
		# 
		# condition: delta_1(source) = a_coefficient (X) b_coefficient (X) (m_coefficient * target)
		def __init__(self, source, target, a_coefficient, b_coefficient, m_coefficient)
			self.source = source
			self.target = target
			self.a_coefficient = a_coefficient
			self.b_coefficient = b_coefficient
			self.m_coefficient = m_coefficient
