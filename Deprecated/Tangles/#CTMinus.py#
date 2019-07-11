from Functions import *

# a class to represent CT^-
# currently, calling CTMinus(tangle)
class CTMinus:
    def __init__(self, tangle):
        self.tangle = tangle

        points = []
        for etangle in tangle.etangles:
            points.append(etangle.left_points())
            points.append(etangle.middle_points())
        points.append(tangle.etangles[-1].right_points())

        self.gens = CTMinus.enumerate_gens(points)

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
