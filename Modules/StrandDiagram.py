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

    def reflect(self):
        return StrandDiagram(
            {orange: tuple(pos[::-1]) for orange, pos in self.orange_strands.items()},
            {orange: sign for orange, sign in self.orange_signs.items()},
            {black: tuple(pos[::-1]) for black, pos in self.black_strands.items()}
        )

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

    def figure_8_case_1b(self, b1: int, b2: int) -> Optional[Dict]:
        a1 = self.black_right_pos(b1)
        a2 = self.black_right_pos(b2)

        for b3 in self.black_strands:
            if not b1 < b3 < b2:
                continue
            a3 = self.black_right_pos(b3)
            if a3 is None:
                return None
            if a3 < a1 or a3 > a2:
                return None

        powers = {}
        for orange in self.orange_strands:
            powers[orange] = 0
            if b1 < self.orange_middle_pos(orange) < b2:
                if self.orange_signs[orange] == 1:
                    if self.orange_left_pos(orange) is not None:  # ???
                        return None
                    if self.orange_right_pos(orange) is not None and \
                            self.orange_right_pos(orange) > a2:
                        powers[orange] += 1
                    if self.orange_right_pos(orange) is not None and \
                            self.orange_right_pos(orange) < a1:
                        powers[orange] += 1
                else:
                    if self.orange_left_pos(orange) is not None:  # ???
                        powers[orange] += 1
                    if self.orange_right_pos(orange) is not None and \
                            self.orange_right_pos(orange) > a2:
                        return None
                    if self.orange_right_pos(orange) is not None and \
                            self.orange_right_pos(orange) < a1:
                        return None

        return powers

    def figure_8_case_2b(self, b1: int, b2: int) -> Optional[Dict]:
        a1 = self.black_left_pos(b1)
        a2 = self.black_left_pos(b2)

        for b3 in self.black_strands:
            if not b1 < b3 < b2:
                continue
            a3 = self.black_left_pos(b3)
            if a3 is None:
                return None
            if a3 < a2 or a3 > a1:
                return None

        powers = {}
        for orange in self.orange_strands:
            powers[orange] = 0
            if b1 < self.orange_middle_pos(orange) < b2:
                if self.orange_signs[orange] == 1:
                    if self.orange_right_pos(orange) is not None:  # ???
                        powers[orange] += 1
                    if self.orange_left_pos(orange) is not None and \
                            self.orange_left_pos(orange) > a1:
                        return None
                    if self.orange_left_pos(orange) is not None and \
                            self.orange_left_pos(orange) < a2:
                        return None
                else:
                    if self.orange_right_pos(orange) is not None:  # ???
                        return None
                    if self.orange_left_pos(orange) is not None and \
                            self.orange_left_pos(orange) > a1:
                        powers[orange] += 1
                    if self.orange_left_pos(orange) is not None and \
                            self.orange_left_pos(orange) < a2:
                        powers[orange] += 1

        return powers

    def figure_8_case_3b(self, b1: int, b2: int) -> Optional[Dict]:
        a1 = self.black_right_pos(b1)
        a2 = self.black_left_pos(b2)

        for b3 in self.black_strands:
            if not b1 < b3 < b2:
                continue
            a3 = self.black_left_pos(b3)
            if a3 is not None:
                if a3 < a2:
                    return None
            else:
                a3 = self.black_right_pos(b3)
                if a3 < a1:
                    return None

        powers = {}
        for orange in self.orange_strands:
            powers[orange] = 0
            if b1 < self.orange_middle_pos(orange) < b2:
                if self.orange_signs[orange] == 1:
                    if self.orange_left_pos(orange) is not None and \
                            self.orange_left_pos(orange) < a2:
                        return None
                    if self.orange_right_pos(orange) is not None and \
                            self.orange_right_pos(orange) < a1:
                        powers[orange] += 1
                else:
                    if self.orange_left_pos(orange) is not None and \
                            self.orange_left_pos(orange) < a2:
                        powers[orange] += 1
                    if self.orange_right_pos(orange) is not None and \
                            self.orange_right_pos(orange) < a1:
                        return None

        return powers

    def figure_8_case_4b(self, b1: int, b2: int) -> Optional[Dict]:
        a1 = self.black_left_pos(b1)
        a2 = self.black_right_pos(b2)

        for b3 in self.black_strands:
            if not b1 < b3 < b2:
                continue
            a3 = self.black_left_pos(b3)
            if a3 is not None:
                if a3 > a1:
                    return None
            else:
                a3 = self.black_right_pos(b3)
                if a3 > a2:
                    return None

        powers = {}
        for orange in self.orange_strands:
            powers[orange] = 0
            if b1 < self.orange_middle_pos(orange) < b2:
                if self.orange_signs[orange] == 1:
                    if self.orange_left_pos(orange) is not None and \
                            self.orange_left_pos(orange) > a1:
                        return None
                    if self.orange_right_pos(orange) is not None and \
                            self.orange_right_pos(orange) > a2:
                        powers[orange] += 1
                else:
                    if self.orange_left_pos(orange) is not None and \
                            self.orange_left_pos(orange) > a1:
                        powers[orange] += 1
                    if self.orange_right_pos(orange) is not None and \
                            self.orange_right_pos(orange) > a2:
                        return None

        return powers

    def figure_8_case_1a(self, b1: int, b2: int) -> Optional[Dict]:
        return self.reflect().figure_8_case_1b(b1, b2)

    def figure_8_case_2a(self, b1: int, b2: int) -> Optional[Dict]:
        return self.reflect().figure_8_case_2b(b1, b2)

    def figure_8_case_3a(self, b1: int, b2: int) -> Optional[Dict]:
        return self.reflect().figure_8_case_3b(b1, b2)

    def figure_8_case_4a(self, b1: int, b2: int) -> Optional[Dict]:
        return self.reflect().figure_8_case_4b(b1, b2)

    # counts the number of orange-black crossings
    def num_orange_black_crossings(self) -> int:
        return sum(self.orange_times_crossed_black(orange, black)
                   for orange in self.orange_strands for black in self.black_strands)

    def black_crosses_black(self, b1: int, b2: int) -> bool:
        return self.black_times_crossed_black(b1, b2) > 0

    def orange_crosses_black(self, orange: int, black: int) -> bool:
        return self.orange_times_crossed_black(orange, black) > 0

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

    def orange_times_crossed_orange(self, o1: int, o2: int) -> int:
        p1 = self.orange_left_pos(o1)
        q1 = self.orange_middle_pos(o1)
        r1 = self.orange_right_pos(o1)
        p2 = self.orange_left_pos(o2)
        q2 = self.orange_middle_pos(o2)
        r2 = self.orange_right_pos(o2)

        return StrandDiagram.times_crossed(p1, q1, r1, p2, q2, r2)

    @staticmethod
    def times_crossed(p1: int, q1: int, r1: int, p2: int, q2: int, r2: int) -> int:
        out = 0
        if None not in {p1,p2,q1,q2}:
            if (p1 < p2) ^ (q1 < q2):
                out += 1
        if None not in {q1,q2,r1,r2}:
            if (q1 < q2) ^ (r1 < r2):
                out += 1

        return out


# more crossing counts to do gradings
    def num_orange_black_pos_crossings(self) -> int:
        return sum(self.orange_times_crossed_black(orange,black)
                   for orange in self.orange_strands if self.orange_signs[orange]==1 for black in self.black_strands)
    def num_orange_black_neg_crossings(self) -> int:
        return sum(self.orange_times_crossed_black(orange,black)
                   for orange in self.orange_strands if self.orange_signs[orange]==-1 for black in self.black_strands)
    def num_orange_orange_pos_crossings(self) -> int:
        return sum(self.orange_times_crossed_orange(o1,o2)
                   for o1 in self.orange_strands if self.orange_signs[o1]==1 for o2 in self.orange_strands if self.orange_signs[o1]==1)
    def num_orange_orange_neg_crossings(self) -> int:
        return sum(self.orange_times_crossed_orange(o1,o2)
                   for o1 in self.orange_strands if self.orange_signs[o1]==-1 for o2 in self.orange_strands if self.orange_signs[o2]==-1)
    def num_black_black_crossings(self) -> int:
        return sum(self.black_times_crossed_black(b1,b2)
                   for b1 in self.black_strands for b2 in self.black_strands) // 2
    def num_negative_orange(self) -> int:
        return sum(1 for o in self.orange_strands if self.orange_signs[o]==-1)

    def nobpc(self) -> int:
        return self.num_orange_black_pos_crossings()
    def nobnc(self) -> int:
        return self.num_orange_black_neg_crossings()
    def noopc(self) -> int:
        return self.num_orange_orange_pos_crossings()
    def noonc(self) -> int:
        return self.num_orange_orange_neg_crossings()
    def nbbc(self) -> int:
        return self.num_black_black_crossings()
    def nno(self) -> int:
        return self.num_negative_orange()



    def maslov(self) -> int:
        left_black = {k:self.black_strands[k] for k in self.black_strands.keys() if self.black_strands[k][2]==None}
        right_black = {k:self.black_strands[k] for k in self.black_strands.keys() if self.black_strands[k][0]==None}
        left_orange = {k:[self.orange_strands[k][0],self.orange_strands[k][1],None]
                       for k in self.orange_strands.keys() if self.orange_strands[k][0] != None}
        right_orange = {k:[None,self.orange_strands[k][1],self.orange_strands[k][2],None]
                        for k in self.orange_strands.keys() if self.orange_strands[k][2] != None}
        left = StrandDiagram(left_orange,self.orange_signs,left_black)
        right = StrandDiagram(right_orange,self.orange_signs,right_black)
        m = right.nbbc()-right.nobpc()+right.noopc()-left.nbbc()+left.nobnc()-left.noonc()-left.nno()
        return m

    def twoalexander(self) -> int:
        twoa = self.nobnc()-self.nobpc()+self.noopc()-self.noonc()-self.nno()
        return twoa



