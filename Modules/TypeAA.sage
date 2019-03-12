from collections import defaultdict
load('Modules/Edge.sage')

class TypeAA:
	# gens - a list of generators
	# edges - a list of Edges
	def __init__(self, gens, edges_out):
		self.gens = gens
		self.edges_out = defaultdict(list,edges_out)
		
	# two cases: other is a TypeDA, or other is a TypeDD
	def tensor(self, other):
		# TODO
		break
		
	def reduce(self):
		# TODO
		break
		