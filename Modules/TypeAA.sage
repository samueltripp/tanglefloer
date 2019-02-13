class TypeAA:
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
		
	# represents the action of some m_{1,i,j} on some pair of generators
	class Edge:
		# source, target - elements of gens
		# l_multipliers - a tuple of i elements of A
		# r_multipliers - a tuple of j elements of B
		# m_cofficient - element of k
		#
		# condition: m_{1,i,j}(source (X) l_multipliers (X) r_multipliers) = m_coefficient * target
		def __init__(self, source, target, l_multipliers, r_multipliers, m_coefficient)
			self.source = source
			self.target = target
			self.l_multipliers = l_multipliers
			self.r_multipliers = r_multipliers
			self.m_coefficient = m_coefficient
