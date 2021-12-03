from __future__ import annotations
from enum import auto, Enum
from typing import Optional, List, Tuple

from SignAlgebra.AMinus import AMinus
from SignAlgebra.Z2PolynomialRing import Z2PolynomialRing, Z2Polynomial


# a tangle
class Tangle:
    # etangles - a list of elementary tangles
    # check - whether or not to check if etangles are compatible with each other
    def __init__(self, etangles):
        for i in range(len(etangles) - 1):
            assert etangles[i].right_sign_sequence() == etangles[i + 1].left_sign_sequence(), \
                "Signs do not match at index {}".format(i)

        self.etangles = tuple(etangles)

        self.height = max(len(etangle.signs) for etangle in etangles)
        self.left_algebra = AMinus(self.left_sign_sequence())
        self.right_algebra = AMinus(self.right_sign_sequence())

    def right_points(self) -> List:
        return self.etangles[-1].right_points()

    # returns the sign sequence corresponding to the left edge of this tangle
    def left_sign_sequence(self):
        return self.etangles[0].left_sign_sequence()

    # returns the sign sequence corresponding to the right edge of this tangle
    def right_sign_sequence(self):
        return self.etangles[-1].right_sign_sequence()

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
        STRAIGHT = auto()

    # etype - one of OVER, UNDER, CUP, CAP, or STRAIGHT
    # signs - a tuple {-1,1}* representing the signs *in the middle* of this tangle
    #         signs[0] is None to enforce 1-indexing
    # position - nat n, where the cup/cap/crossing is between strand indices n and n+1
    #           None if etype is STRAIGHT
    # over/under represents what the bottom strand does; under on the left, over on the right
    # CONVENTION: strand indices are 1-indexed, starting from the bottom
    def __init__(self, etype: ETangle.Type, signs, position=None):

        for sign in signs:
            assert sign in (-1, 1), "{} is not a valid sign.".format(sign)
        if position is not None:
            assert 1 <= position < len(signs)+1, "{} is not a valid position.".format(position)

        if etype in (ETangle.Type.CUP, ETangle.Type.CAP):
            assert signs[position - 1] == -signs[position], "Signs are not compatible with {}.".format(etype)

        self.etype = etype
        self.signs = (None,) + signs
        self.position = position

        super().__init__((self,))

        self.ring = \
            Z2PolynomialRing(['U' + str(i)
                              for i in range(1, len(self.middle_points()) + len(self.left_algebra.positives) - 1)])

        self.left_scalar_action = self.build_left_scalar_action()
        self.right_scalar_action = self.build_right_scalar_action()

    def strand_index_to_variable_name(self, strand_index: int) -> str:
        if self.etype in (ETangle.Type.CUP, ETangle.Type.CAP) and strand_index > self.position:
            return 'U' + str(strand_index-1)
        return 'U' + str(strand_index)

    def left_strand_position_to_index(self, pos):
        if self.etype in (ETangle.Type.CAP, ETangle.Type.OVER, ETangle.Type.STRAIGHT):
            return pos
        elif self.etype == ETangle.Type.CUP:
            if pos >= self.position + 1:
                return pos - 1
            else:
                return pos
        elif self.etype == ETangle.Type.UNDER:
            if pos == self.position:
                return self.position + 1
            elif pos == self.position + 1:
                return self.position
            else:
                return pos
        else:
            raise Exception('unknown ETangle.Type')

    def right_strand_position_to_index(self, pos):
        if self.etype in (ETangle.Type.CUP, ETangle.Type.UNDER, ETangle.Type.STRAIGHT):
            return pos
        elif self.etype == ETangle.Type.CAP:
            if pos >= self.position + 1:
                return pos - 1
            else:
                return pos
        elif self.etype == ETangle.Type.OVER:
            if pos == self.position:
                return self.position + 1
            elif pos == self.position + 1:
                return self.position
            else:
                return pos
        else:
            raise Exception('unknown ETangle.Type')

    # turns the given strand index into a variable in F[U1, U2, ...]
    def strand_index_to_variable(self, strand_index: int) -> Z2Polynomial:
        return self.ring[self.strand_index_to_variable_name(strand_index)]

    def left_algebra_strand_index_to_variable(self, strand_index: int) -> Z2Polynomial:
        return self.left_algebra.ring['U' +
                                      str(self.left_algebra.positives.index(self.left_y_pos(strand_index) + 1 / 2))]

    def build_left_scalar_action(self):
        return Z2PolynomialRing.Map(self.left_algebra.ring, self.ring,
                                    {'U' + str(i): 'U' + str(len(self.middle_points()) - 1 + i)
                                     for i, p in enumerate(self.left_algebra.positives)
                                        if p is not None})

    def build_right_scalar_action(self):
        return Z2PolynomialRing.Map(self.right_algebra.ring, self.ring,
                                    {'U' + str(i): self.strand_index_to_variable_name(
                                        self.right_strand_position_to_index(p))
                                     for i, p in enumerate(self.right_algebra.positives)
                                        if p is not None and self.right_strand_position_to_index(p) is not None})

    def strand_index_to_left_sign(self, strand_index) -> Optional[int]:
        if self.etype == ETangle.Type.CUP and strand_index in (self.position, self.position + 1):
            return None
        else:
            return self.signs[strand_index]

    # returns the sign sequence corresponding to the left edge of this tangle
    def left_sign_sequence(self) -> Tuple:
        if self.etype == ETangle.Type.CUP:
            return self.signs[:self.position] + \
                   self.signs[self.position + 2:]
        elif self.etype == ETangle.Type.UNDER:
            return self.signs[:self.position] + \
                   (self.signs[self.position+1], self.signs[self.position]) + \
                   self.signs[self.position+2:]
        else:
            return self.signs

    # returns the sign sequence corresponding to the right edge of this tangle
    def right_sign_sequence(self) -> Tuple:
        if self.etype == ETangle.Type.CAP:
            return self.signs[:self.position] + self.signs[self.position + 2:]
        elif self.etype == ETangle.Type.OVER:
            return self.signs[:self.position] + \
                   (self.signs[self.position+1], self.signs[self.position]) + \
                   self.signs[self.position+2:]
        else:
            return self.signs

    # given a strand index, returns the y-position of that strand on the left
    def left_y_pos(self, strand_index: int) -> Optional[float]:
        if self.etype in (ETangle.Type.OVER, ETangle.Type.CAP, ETangle.Type.STRAIGHT):
            return strand_index - 1/2
        elif self.etype == ETangle.Type.UNDER:
            if strand_index == self.position:
                return self.position + 1/2
            elif strand_index == self.position + 1:
                return self.position - 1/2
            else:
                return strand_index - 1/2
        elif self.etype == ETangle.Type.CUP:
            if strand_index < self.position:
                return strand_index - 1/2
            elif strand_index > self.position + 1:
                return strand_index - 5/2
            else:
                return None
        else:
            raise Exception('unknown ETangle.Type')

    # given a strand index, returns the y-position of that strand in the middle
    def middle_y_pos(self, strand_index: int) -> float:
        if self.etype in (ETangle.Type.OVER, ETangle.Type.UNDER, ETangle.Type.STRAIGHT):
            return strand_index - 1/2
        else:
            if strand_index in (self.position, self.position + 1):
                return self.position
            else:
                return strand_index - 1/2

    # given a strand index, returns the y-position of that strand on the right
    def right_y_pos(self, strand_index: int) -> Optional[float]:
        if self.etype in (ETangle.Type.UNDER, ETangle.Type.CUP, ETangle.Type.STRAIGHT):
            return strand_index - 1 / 2
        elif self.etype == ETangle.Type.OVER:
            if strand_index == self.position:
                return self.position + 1 / 2
            elif strand_index == self.position + 1:
                return self.position - 1 / 2
            else:
                return strand_index - 1 / 2
        elif self.etype == ETangle.Type.CAP:
            if strand_index < self.position:
                return strand_index - 1 / 2
            elif strand_index > self.position + 1:
                return strand_index - 5 / 2
            else:
                return None
        else:
            raise Exception('unknown ETangle.Type')

    # does this strand stay straight on the left side of this tangle?
    def left_strand_straight(self, strand_index: int) -> bool:
        return self.left_y_pos(strand_index) == self.middle_y_pos(strand_index)

    # does this strand stay straight on the right side of this tangle?
    def right_strand_straight(self, strand_index: int) -> bool:
        return self.right_y_pos(strand_index) == self.middle_y_pos(strand_index)

    # returns the list of points corresponding to the left side of this tangle
    def left_points(self) -> List:
        if self.etype == ETangle.Type.CUP:
            return list(range(len(self.signs) - 2))
        return list(range(len(self.signs)))

    # returns the list of points corresponding to the middle of this tangle
    def middle_points(self) -> List:
        if self.etype in (ETangle.Type.CUP, ETangle.Type.CAP):
            return list(range(self.position)) + list(range(self.position + 1, len(self.signs)))
        return list(range(len(self.signs)))

    # returns the list of points corresponding to the right side of this tangle
    def right_points(self) -> List:
        if self.etype == ETangle.Type.CAP:
            return list(range(len(self.signs) - 2))
        return list(range(len(self.signs)))

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
