import numpy as np
from collections import defaultdict
load('Modules/TypeDD.sage')
load('Modules/TypeAA.sage')
load('Modules/Edge.sage')

class TypeDA:
	# gens - a list of generators
	# edges_out - a dictionary of edges keyed by source generator
	# dictionary values are lists of edges (each edge is a 5-tuple (source, target, a, b, m))

	def __init__(self, gens, edges_out):
		self.gens = gens
		self.edges_out = defaultdict(list,edges_out)
	
	# only to be paired with 'other' with class TypeDD. The code can be modified to also pair with type DA. In this case, the label of b_coefficient vs. b_coefficient needs to be modified, and the output of the program needs to be changed to typeDA
	def tensor(self, other):
		#the following code makes a matrix, indexed by the DA and DD gens resp., which contains idempotent data for the box tensor product generators. Entry (i,j) is zero if the ith DA generator does not pair with the jth DD generator.
		MNgens = np.zeros((len(self.gens),len(other.gens)),dtype=object)
		listkey = 0
		for i in range(len(self.gens)):
			for j in range(len(other.gens)):
				if self.gens[i][1] == other.gens[j][0]: 
					MNgens[i][j] = [self.gens[i][0],other.gens[j][1],listKey]
					listKey += 1 
					
		#the following code makes edges for the tensor product.
		MNedgeDict = defaultdict(list)  

		for i in range(len(self.gens)):
			for e in self.edges_out[i]:
				for j in range(len(other.gens)):
					if MNgens[i][j] != 0:
						if e.b_coefficient == []: #this is the case where the DA edge has a trivial right multiplier, in which case we pair with an identity in the DD module
							MNedgeDict[MNgens[i][j][2]].append(Edge(MNgens[i][j][2], MNgens[e.target][j][2], e.a_coefficient, [], e.m_coefficient))
						else:  #this is the case where the multipliers on the right of the DA edge are nontrivial, and therefore need to pair with the left coefficient of some DD module edge. 
							#TODO: this code needs to be modified if multiple e.b_coefficient are allowed. We will need to search through a path of DD module edges to pair with it in this case.
							for ee in other.edges_out[j]: 
								if e.b_coefficient ==  ee.a_coefficient:
									MNedgeDict[MNgens[i][j][2]].append(Edge(MNgens[i][j][2],MNgens[e.target][ee.target][2], e.a_coefficient,ee.b_coefficient, e.m_coefficient*ee.m_coefficient))

		#this is the case where a DD edge pairs with an identity edge in the DA module. 
		for j in range(len(other.gens)):
			for ee in other.edges_out[j]:
				if ee.a_coefficient == []:
					for i in range(len(self.gens)):
						if MNgens[i][j] != 0: 
							MNedgeDict[MNgens[i][j][2]].append(Edge(MNgens[i][j][2],MNgens[i][ee.target][2], [],ee.b_coefficient, ee.m_coefficient))

		#this makes a list of the tensor product generators
		outGenList = []
		for i in range(len(self.gens)):
			for j in range(len(other.gens)):
				if MNgens[i][j] != 0:
					outGenList.append([MNgens[i][j][0],MNgens[i][j][1]])

		#the box tensor of a DA with a DD module results in a type DD module
		return TypeDD(outGenList,MNedgeDict)				

		
	def reduce(self):
		# TODO
		return 0
		