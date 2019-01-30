from enum import Enum
import numpy


# a tangle
class Tangle:
    # etangles - a list of elementary tangles
    # check - whether or not to check if etangles are compatible with each other
    def __init__(self, etangles, check=True):
        if check:
            for i in range(len(etangles) - 1):
                assert etangles[i].right_signs()== etangles[i + 1].left_signs(), \
                    "Signs do not match at index {}".format(i)

        self.etangles = etangles

    def left_signs(self):
        return self.etangles[0].left_signs()

    def right_signs(self):
        return self.etangles[-1].right_signs()

    def __add__(self, other):
        assert self.right_signs() == other.left_signs(), "Signs do not match."
        return Tangle(self.etangles + other.etangles, False)

    def __eq__(self, other):
        if isinstance(other, Tangle):
            return self.etangles == other.etangles
        else:
            return False


# an elementary tangle
class ETangle(Tangle):
    Type = Enum('Type', 'OVER UNDER CUP CAP')

    # etype - one of OVER, UNDER, CUP, or CAP
    # signs - a tuple {-1,1}* representing the signs of the left edge of the tangle,
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

        self.etangles = [self]

    def left_signs(self):
        if self.etype == ETangle.Type.CUP:
            return self.signs[:self.position - 1] + self.signs[self.position + 1:]
        else:
            return self.signs

    def right_signs(self):
        if self.etype == ETangle.Type.CAP:
            return self.signs[:self.position - 1] + self.signs[self.position + 1:]
        elif self.etype in (ETangle.Type.OVER, ETangle.Type.UNDER):
            return self.signs[:self.position - 1] + (
            self.signs[self.position], self.signs[self.position - 1]) + self.signs[self.position + 1:]
        else:
            return self.signs

    def __eq__(self, other):
        if isinstance(other, Tangle) and len(other.etangles)==1:
            other_etangle = other.etangles[0]
            return self.etype == other_etangle.etype \
                and self.signs == other_etangle.signs \
                and self.position == other_etangle.position
        else:
            return False