from enum import Enum
import numpy


# an elementary tangle
class ETangle:
    Type = Enum('Type', 'OVER UNDER CUP CAP')

    # etype - one of OVER, UNDER, CUP, or CAP
    # signs - a tuple {-1,1}* representing the signs on the left edge of the tangle,
    #         unless etype is CUP, in which case it represents the signs on the right edge
    # position - nat n, where the cup/cap/crossing is between places n-1 and n
    def __init__(self, etype, signs, position):
        assert etype in list(ETangle.Type), "{} is not a valid type.".format(etype)
        for sign in signs:
            assert sign in (-1, 1), "{} is not a valid sign.".format(sign)
        assert 0 < position < len(signs), "{} is not a valid position.".format(position)

        if etype in (ETangle.Type.CUP, ETangle.Type.CAP):
            assert signs[position - 1] == -signs[position], "Signs are not compatible with {}.".format(etype)

        self.etype = etype
        self.signs = signs
        self.position = position

    # returns the signs on the left edge of the tangle
    def left_signs(self):
        if self.etype == ETangle.Type.CUP:
            return self.signs[:self.position - 1] + self.signs[self.position + 1:]
        else:
            return self.signs

    # returns the signs on the right edge of the tangle
    def right_signs(self):
        if self.etype == ETangle.Type.CAP:
            return self.signs[:self.position - 1] + self.signs[self.position + 1:]
        elif self.etype in (ETangle.Type.OVER, ETangle.Type.UNDER):
            return self.signs[:self.position - 1] + (
            self.signs[self.position], self.signs[self.position - 1]) + self.signs[self.position + 1:]
        else:
            return self.signs

    _ascii_cup = numpy.array([list(l) for l in [' | ', '` ,', '- -', '- -']])
    _ascii_cap = numpy.array([list(l) for l in ['- -', '- -', '\' .', ' | ']])
    _ascii_over = numpy.array([list(l) for l in ['- -', '- -', '- -', '- -', '\' .', ' / ', '` ,', '- -']])
    _ascii_under = numpy.array([list(l) for l in ['- -', '\' .', ' \ ', '` ,', '- -', '- -', '- -', '- -']])

    def ascii_array(self, max_strands):
        a = numpy.array([[' ' for y in range(2 * max_strands)] for x in range(8)])
        for i in range(self.position - 1):
            for x in range(8):
                a[x, 2 * i] = '-'
        if self.etype in (ETangle.Type.CUP, ETangle.Type.CAP):
            if self.etype == ETangle.Type.CUP:
                a[4:, 2 * self.position - 2:2 * self.position + 1] = ETangle._ascii_cup
                for i in range(self.position + 1, len(self.signs)):
                    a[0, 2 * i - 4] = '\''
                    for x in range(1, 4):
                        a[x, 2 * i + x - 4] = '/'
                    a[4, 2 * i] = ','
                    for x in range(5, 8):
                        a[x, 2 * i] = '-'
            else:
                a[:4, 2 * self.position - 2:2 * self.position + 1] = ETangle._ascii_cap
                for i in range(self.position + 1, len(self.signs)):
                    for x in range(0, 3):
                        a[x, 2 * i] = '-'
                    a[3, 2 * i] = '.'
                    for x in range(4, 7):
                        a[x, 2 * i - x + 3] = '\\'
                    a[7, 2 * i - 4] = '`'
        else:
            for i in range(self.position + 1, len(self.signs)):
                for x in range(8):
                    a[x, 2 * i] = '-'
            if self.etype == ETangle.Type.OVER:
                a[:, 2 * self.position - 2:2 * self.position + 1] = ETangle._ascii_over
            else:
                a[:, 2 * self.position - 2:2 * self.position + 1] = ETangle._ascii_under
        return a


# a - a numpy array of characters, indexed [x,y]
def draw_ascii_array(a):
    return ''.join([''.join(a[:, y]) + '\n' for y in range(len(a[0]) - 1, -1, -1)])


# a tangle
class Tangle:
    # etangles - a sequence of elementary tangles
    def __init__(self, etangles):
        for i in range(len(etangles) - 1):
            assert etangles[i].right_signs() == etangles[i + 1].left_signs(), "Signs do not match at index {}".format(i)

        self.etangles = etangles

    # draws an ascii diagram of the tangle
    def draw_ascii(self):
        max_strands = max([len(etangle.signs) for etangle in self.etangles])
        return draw_ascii_array(numpy.concatenate([etangle.ascii_array(max_strands) for etangle in self.etangles]))