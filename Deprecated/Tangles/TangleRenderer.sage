import numpy
# import svgwrite
load('Tangle.sage')
load('Functions.sage')


# utility class for drawing tangles
class TangleRenderer:
	_CUP = numpy.array([list(l) for l in [' | ', '` ,', '- -', '- -']])
	_CAP = numpy.array([list(l) for l in ['- -', '- -', '\' .', ' | ']])
	_OVER = numpy.array([list(l) for l in ['- -', '- -', '- -', '- -', '\' .', ' / ', '` ,', '- -']])
	_UNDER = numpy.array([list(l) for l in ['- -', '\' .', ' \ ', '` ,', '- -', '- -', '- -', '- -']])

	# turn the given etangle into a numpy array of characters with the given strand height
	@staticmethod
	def ascii_array(etangle, max_strands):
		a = numpy.array([[' ' for y in range(2 * max_strands)] for x in range(8)])
		for i in range(etangle.position - 1):
			for x in range(8):
				a[x, 2 * i] = '-'
			a[3, 2 * i] = '>' if etangle.signs[i]==1 else '<'
		if etangle.etype in (ETangle.Type.CUP, ETangle.Type.CAP):
			if etangle.etype == ETangle.Type.CUP:
				a[4:, 2 * etangle.position - 2:2 * etangle.position + 1] = TangleRenderer._CUP
				for i in range(etangle.position + 1, len(etangle.signs)):
					a[0, 2 * i - 4] = '\''
					for x in range(1, 4):
						a[x, 2 * i + x - 4] = '/'
					a[4, 2 * i] = ','
					for x in range(5, 8):
						a[x, 2 * i] = '-'
			else:
				a[:4, 2 * etangle.position - 2:2 * etangle.position + 1] = TangleRenderer._CAP
				for i in range(etangle.position + 1, len(etangle.signs)):
					for x in range(0, 3):
						a[x, 2 * i] = '-'
					a[3, 2 * i] = '.'
					for x in range(4, 7):
						a[x, 2 * i - x + 3] = '\\'
					a[7, 2 * i - 4] = '`'
		else:
			for i in range(etangle.position + 1, len(etangle.signs)):
				for x in range(8):
					a[x, 2 * i] = '-'
				a[3, 2 * i] = '>' if etangle.signs[i]==1 else '<'
			if etangle.etype == ETangle.Type.OVER:
				a[:, 2 * etangle.position - 2:2 * etangle.position + 1] = TangleRenderer._OVER
			else:
				a[:, 2 * etangle.position - 2:2 * etangle.position + 1] = TangleRenderer._UNDER
		return a

	# a - a numpy array of characters, indexed [x,y]
	@staticmethod
	def ascii_array_to_string(a):
		return ''.join([''.join(a[:, y]) + '\n' for y in range(len(a[0]) - 1, -1, -1)])

	# create an ascii string from the given tangle
	@staticmethod
	def ascii(tangle):
		max_strands = max([len(etangle.signs) for etangle in tangle.etangles])
		return TangleRenderer.ascii_array_to_string(
			numpy.concatenate([TangleRenderer.ascii_array(etangle, max_strands) for etangle in tangle.etangles]))

#     # create an SVG diagram from the given elementary tangle
#     @staticmethod
#     def svg_helper(etangle):
#         paths = []
#
#         if etangle.etype in (ETangle.Type.OVER, ETangle.Type.UNDER):
#             for i in range(etangle.position - 1):
#                 paths.append(line(c(0, i), (1, i)))
#             for i in range(etangle.position + 1, len(etangle.signs)):
#                 paths.append(line(c(0, i), c(1, i)))
#             if etangle.etype == ETangle.Type.OVER:
#                 paths.append(line(c(0, etangle.position-1), c(0.5, etangle.position-1)))
#                 paths.append(line(c(0, etangle.position), c(0.5, etangle.position)))
#                 paths.append(
#                     cubic_bezier(c(0.5, etangle.position - 1), c(0.75, etangle.position - 1), c(0.75, etangle.position),
#                                 c(1, etangle.position)))
#                 paths.append(cubic_bezier(c(0.5, etangle.position), c(0.75, etangle.position), c(0.75, etangle.position - 1),
#                                 c(1, etangle.position - 1), under=True))
#             elif etangle.etype == ETangle.Type.UNDER:
#                 paths.append(
#                     cubic_bezier(c(0, etangle.position), c(0.25, etangle.position), c(0.25, etangle.position - 1),
#                                 c(0.5, etangle.position - 1)))
#                 paths.append(cubic_bezier(c(0, etangle.position - 1), c(0.25, etangle.position - 1), c(0.25, etangle.position),
#                                 c(0.5, etangle.position), under=True))
#                 paths.append(line(c(0.5, etangle.position-1), c(1, etangle.position-1)))
#                 paths.append(line(c(0.5, etangle.position-1), c(1, etangle.position-1)))
#                 paths.append(line(c(0.5, etangle.position), c(1, etangle.position)))
#         elif etangle.etype == ETangle.Type.CUP:
#             for i in range(etangle.position - 1):
#                 paths.append(line(c(0, i), c(1, i)))
#             for i in range(etangle.position + 1, len(etangle.signs)):
#                 paths.append(cubic_bezier(c(0, i), c(0.25, i), c(0.25, i + 2), c(0.5, i + 2)))
#                 paths.append(line(c(0.5, i), c(1, i)))
#             paths.append(semicircle(c(1, etangle.position - 1), 0.5, 1))
#         elif etangle.etype == ETangle.Type.CAP:
#             for i in range(etangle.position - 1):
#                 paths.append(line(c(0, i), c(1, i)))
#             for i in range(etangle.position + 1, len(etangle.signs)):
#                 paths.append(line(c(0, i), c(0.5, i)))
#                 paths.append(cubic_bezier(c(0.5, i), c(0.75, i), c(0.75, i - 2), c(1, i - 2)))
#             paths.append(semicircle(c(0, etangle.position-1), 0.5, 0))
#
#         return paths
#
#     @staticmethod
#     def svg(filename, tangle, generator=None):
#         paths = []
#         for i, etangle in enumerate(tangle.etangles):
#             for path in TangleRenderer.svg_helper(etangle):
#                 path.translate(i,0)
#                 paths.append(path)
#         if generator:
#             for i, partial_bijection in enumerate(generator):
#                 for x,y in partial_bijection.items():
#                     paths.append(cubic_bezier(c(i/2,x-.5), c(i/2+0.25,x-.5), c(i/2+0.25, y-.5), c(i/2+0.5,y-.5), 'black'))
#         dwg = svgwrite.Drawing(filename=filename, debug=True)
#         for path in paths:
#             dwg.add(path)
#         dwg.viewbox(0,-max([len(etangle.signs) for etangle in tangle.etangles]),len(tangle.etangles)+0.5,max([len(etangle.signs) for etangle in tangle.etangles])+1)
#         dwg.save(pretty=True)
#
#
# def cubic_bezier(x1, x2, x3, x4, color='orange', under = False):
#     if under:
#         return svgwrite.path.Path('M {} {} C {} {} {} {} {} {}'.format(*x1, *x2, *x3, *x4), stroke = color, stroke_width = '0.05', fill='none', stroke_dasharray = '0.4')
#     else:
#         return svgwrite.path.Path('M {} {} C {} {} {} {} {} {}'.format(*x1, *x2, *x3, *x4), stroke = color, stroke_width = '0.05', fill='none')
#
#
# def line(x1, x2, color='orange'):
#     return svgwrite.path.Path('M {} {} L {} {}'.format(*x1, *x2), stroke = color, stroke_width = '0.05', fill='none')
#
#
# def semicircle(s, r, direction, color='orange'):
#     return svgwrite.path.Path('M {} {} A {} {} 0 1 {} {} {}'.format(*s, r, r, direction, s[0], s[1]-2*r), stroke = color, stroke_width = '0.05', fill='none')
