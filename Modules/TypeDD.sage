import numpy as np

class TypeDD:
	# gens - a list of generators
	# edges - a list of Edges
	def __init__(self, gens, edges):
		self.gens = gens
		self.edges = dictEdges(edges)
	
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