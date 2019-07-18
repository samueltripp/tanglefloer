from __future__ import annotations
from enum import *
from SignAlgebra.Z2PolynomialRing import *
from SignAlgebra.AMinus import *
import numpy


# a tangle
class Tangle:
    # etangles - a list of elementary tangles
    # check - whether or not to check if etangles are compatible with each other
    def __init__(self, etangles):
        for i in range(len(etangles) - 1):
            assert etangles[i].right_signs() == etangles[i + 1].left_signs(), \
                "Signs do not match at index {}".format(i)

        self.etangles = tuple(etangles)

        height = max(len(etangle.signs) for etangle in etangles)
        self.polyring = Z2PolynomialRing(['U%s' % p for p in range(1, height + 1)])
        self.left_algebra = AMinus(self.left_signs(), self.polyring)
        self.right_algebra = AMinus(self.right_signs(), self.polyring)

    # returns the sign sequence corresponding to the left edge of this tangle
    def left_signs(self):
        return self.etangles[0].left_signs()

    # returns the sign sequence corresponding to the right edge of this tangle
    def right_signs(self):
        return self.etangles[-1].right_signs()

    # add two tangles
    def __add__(self, other):
        return Tangle(self.etangles + other.etangles)

    # right addition is used for things like sum()
    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return other + self

    def __eq__(self, other: Tangle):
        return self.etangles == other.etangles

    def __str__(self):
        from Tangles.TangleRenderer import TangleRenderer
        return TangleRenderer.ascii(self)


# an elementary tangle
class ETangle(Tangle):
    class Type(Enum):
        OVER = auto()
        UNDER = auto()
        CUP = auto()
        CAP = auto()

    # etype - one of OVER, UNDER, CUP, or CAP
    # signs - a tuple {-1,1}* representing the signs of the left edge of the tangle,
    #         unless etype is CUP, in which case it represents the signs on the right edge
    # position - nat n, where the cup/cap/crossing is between places n-1 and n
    # over/under represents what the bottom strand does; under on the left, over on the right
    def __init__(self, etype: ETangle.Type, signs, position):
        for sign in signs:
            assert sign in (-1, 1), "{} is not a valid sign.".format(sign)
        assert 0 < position < len(signs), "{} is not a valid position.".format(position)

        if etype in (ETangle.Type.CUP, ETangle.Type.CAP):
            assert signs[position - 1] == -signs[position], "Signs are not compatible with {}.".format(etype)

        self.etype = etype
        self.signs = signs
        self.position = position

        super().__init__((self,))

    # returns the sign sequence corresponding to the left edge of this tangle
    def left_signs(self):
        if self.etype == ETangle.Type.CUP:
            return self.signs[:self.position - 1] + self.signs[self.position + 1:]
        else:
            return self.signs

    # returns the sign sequence corresponding to the right edge of this tangle
    def right_signs(self):
        if self.etype == ETangle.Type.CAP:
            return self.signs[:self.position - 1] + self.signs[self.position + 1:]
        elif self.etype in (ETangle.Type.OVER, ETangle.Type.UNDER):
            return self.signs[:self.position - 1] + \
                   (self.signs[self.position], self.signs[self.position - 1]) + \
                   self.signs[self.position + 1:]
        else:
            return self.signs

    # returns the set of points corresponding to the left side of this tangle
    def left_points(self):
        if self.etype == ETangle.Type.CUP:
            return list(range(len(self.signs) - 1))
        return list(range(len(self.signs) + 1))

    # returns the set of points corresponding to the middle of this tangle
    def middle_points(self):
        if self.etype in (ETangle.Type.CUP, ETangle.Type.CAP):
            return list(range(self.position)) + list(range(self.position + 1, len(self.signs) + 1))
        return list(range(len(self.signs) + 1))

    # returns the set of points corresponding to the right side of this tangle
    def right_points(self):
        if self.etype == ETangle.Type.CAP:
            return list(range(len(self.signs) - 1))
        return list(range(len(self.signs) + 1))

    # returns the mapping {left_point: middle_point}
    def left_to_middle(self) -> Dict:
        num_left_points = len(self.left_points())
        if self.etype == ETangle.Type.OVER:
            out = {p: p for p in range(num_left_points)}
            return out
        elif self.etype == ETangle.Type.UNDER:
            out = {p: p for p in range(num_left_points)}
            out[self.position - 1] = self.position
            out[self.position] = self.position - 1
            return out
        elif self.etype == ETangle.Type.CUP:
            out = {}
            for p in range(self.position - 1):
                out[p] = p
            for p in range(self.position, num_left_points):
                out[p] = p + 2
            return out
        elif self.etype == ETangle.Type.CAP:
            out = {p: p for p in range(num_left_points)}
            del out[self.position - 1]
            del out[self.position]
            return out

    # returns the mapping {right_point: middle_point}
    def right_to_middle(self) -> Dict:
        num_right_points = len(self.right_points())
        if self.etype == ETangle.Type.UNDER:
            out = {p: p for p in range(num_right_points)}
            return out
        elif self.etype == ETangle.Type.OVER:
            out = {p: p for p in range(num_right_points)}
            out[self.position - 1] = self.position
            out[self.position] = self.position - 1
            return out
        elif self.etype == ETangle.Type.CAP:
            out = {}
            for p in range(self.position - 1):
                out[p] = p
            for p in range(self.position, num_right_points):
                out[p] = p + 2
            return out
        elif self.etype == ETangle.Type.CUP:
            out = {p: p for p in range(num_right_points)}
            del out[self.position - 1]
            del out[self.position]
            return out

    def __repr__(self):
        return str((self.etype, self.signs, self.position))

    def __eq__(self, other):
        if isinstance(other, Tangle) and len(other.etangles) == 1:
            other_etangle = other.etangles[0]
            return self.etype == other_etangle.etype \
                and self.signs == other_etangle.signs \
                and self.position == other_etangle.position
        else:
            return False

    def __hash__(self):
        return hash(self.etype) + hash(self.signs) + hash(self.position)
