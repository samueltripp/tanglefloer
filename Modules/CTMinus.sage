load("Tangles/Functions.sage")

# a class to represent CT^-
# this should contain all the necessary methods for turning a tangle into a type DA structure
class CTMinus:
	def __init__(self, tangle):
		self.tangle = tangle

	def toTypeDA(self):
		return TypeDA(self.generators(), self.edges())
		
	# returns the generators of CT^-(T)
	def generators(self):
		points = []
		for etangle in self.tangle.etangles:
			points.append(etangle.left_points())
			points.append(etangle.middle_points())
		points.append(self.tangle.etangles[-1].right_points())

		return CTMinus.enumerate_gens(points)

	# points - a list of sets of points
	@staticmethod
	def enumerate_gens(points):
		def reverse_injection(d):
			return {v: k for k, v in d.items()}

		l = len(points)
		sequences = []
		if l < 2:
			return sequences
		elif l == 2:
			return [[pb] for pb in partial_bijections(points[0],points[1])]
		else:
			mid = l // 2
			for pb in partial_bijections(points[mid],points[mid + 1]):
				r, s = list(reversed(points[:mid + 1])), points[mid + 1:]

				r[0] = list(set(r[0]).difference(pb.keys()))
				s[0] = list(set(s[0]).difference(pb.values()))

				t = list(CTMinus.enumerate_gens_helper(s))
				sequences.extend([[reverse_injection(inj) for inj in reversed(s1)] \
					+ [pb] + s2 for s1 in CTMinus.enumerate_gens_helper(r) for s2 in t])
			return sequences

	@staticmethod
	def enumerate_gens_helper(points):
		if len(points) == 2:
			for inj in injections(points[0],points[1]):
				yield [inj]
		elif len(points) > 2:
			for inj in injections(points[0],points[1]):
				coker = set(points[1]).difference(inj.values())
				for sequence in CTMinus.enumerate_gens_helper([list(coker)]+points[2:]):
					yield [inj] + sequence

	# should return the list of graph edges corresponding to the maps delta_1^j
	def edges(self):
		# TODO
		# something involving d_plus, d_minus, d_mixed
		return
		
	def d_plus(self):
		# TODO
		return

	def d_minus(self):
		# TODO
		return

	def d_mixed(self):
		# TODO
		return