import numpy
from Tangles import *


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
