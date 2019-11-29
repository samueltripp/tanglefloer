from __future__ import annotations
from typing import Dict, Tuple, Optional
import itertools
from SignAlgebra.Z2PolynomialRing import Z2Polynomial


class StrandDiagram:
    # orange_strands: {strand_index: (left_pos, middle_pos, right_pos)}
    # orange_signs: {strand_index: int}
    # black_strands {strand_index: (left_pos, middle_pos, right_pos)}
    def __init__(self, orange_strands: Dict, orange_signs: Dict, black_strands: Dict):
        self.orange_strands = orange_strands
        self.orange_signs = orange_signs
        self.black_strands = black_strands

    def orange_left_pos(self, orange_index: int):
        return self.orange_strands[orange_index][0]

    def orange_middle_pos(self, orange_index: int):
        return self.orange_strands[orange_index][1]

    def orange_right_pos(self, orange_index: int):
        return self.orange_strands[orange_index][2]

    def black_left_pos(self, orange_index: int):
        return self.black_strands[orange_index][0]

    def black_middle_pos(self, orange_index: int):
        return self.black_strands[orange_index][1]

    def black_right_pos(self, orange_index: int):
        return self.black_strands[orange_index][2]

    # simplify this diagram using the relations from Figure 6 in "An introduction..."
    # output dictionary is {orange_strand_index: power of corresponding variable}
    def figure_6_relations(self) -> Optional[Dict]:
        powers = {}
        for b1, b2 in itertools.combinations(self.black_strands.keys(), 2):
            if self.black_double_crosses_black(b1, b2):
                return None

        for orange in self.orange_strands:
            powers[orange] = 0
            for black in self.black_strands:
                if self.orange_double_crosses_black(orange, black):
                    if self.orange_signs[orange] == 1:
                        powers[orange] += 1
                    else:
                        return None

        return powers

    def figure_7_relations(self) -> Optional[Dict]:
        powers = {}
        for b1, b2 in itertools.combinations(self.black_strands.keys(), 2):
            if self.black_double_crosses_black(b1, b2):
                return None

        for orange in self.orange_strands:
            powers[orange] = 0
            for black in self.black_strands:
                if self.orange_double_crosses_black(orange, black):
                    if self.orange_signs[orange] == -1:
                        powers[orange] += 1
                    else:
                        return None

        return powers

    def figure_8_relations(self) -> Optional[Tuple[Dict, StrandDiagram]]:
        pass  # TODO

    def black_double_crosses_black(self, b1: int, b2: int) -> bool:
        return self.black_times_crossed_black(b1, b2) > 1

    def orange_double_crosses_black(self, orange: int, black: int) -> bool:
        return self.orange_times_crossed_black(orange, black) > 1

    def black_times_crossed_black(self, b1: int, b2: int) -> int:
        p1 = self.black_left_pos(b1)
        q1 = self.black_middle_pos(b1)
        r1 = self.black_right_pos(b1)
        p2 = self.black_left_pos(b2)
        q2 = self.black_middle_pos(b2)
        r2 = self.black_right_pos(b2)

        return StrandDiagram.times_crossed(p1, q1, r1, p2, q2, r2)

    def orange_times_crossed_black(self, orange: int, black: int) -> int:
        p1 = self.orange_left_pos(orange)
        q1 = self.orange_middle_pos(orange)
        r1 = self.orange_right_pos(orange)
        p2 = self.black_left_pos(black)
        q2 = self.black_middle_pos(black)
        r2 = self.black_right_pos(black)

        return StrandDiagram.times_crossed(p1, q1, r1, p2, q2, r2)

    @staticmethod
    def times_crossed(p1: int, q1: int, r1: int, p2: int, q2: int, r2: int) -> int:
        out = 0
        if (p1 < p2) ^ (p1 < q2):
            out += 1
        if (p1 < q2) ^ (p1 < r2):
            out += 1
        if (p1 < r2) ^ (q1 < r2):
            out += 1
        if (q1 < r2) ^ (r1 < r2):
            out += 1

        return out
