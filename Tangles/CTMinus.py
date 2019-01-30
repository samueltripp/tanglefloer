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
        sequences = []
        if len(points) < 2:
            return sequences
        elif len(points) == 2:
            return [[pb] for pb in partial_bijections(points[0],points[1])]
        else:
            for pb in partial_bijections(points[0],points[1]):
                coker = set(points[1]).difference(pb.values())
                sequences.extend([[pb]+sequence for sequence in CTMinus.enumerate_gens_helper([list(coker)]+points[2:])])
        return sequences

    @staticmethod
    def enumerate_gens_helper(points):
        sequences = []
        if len(points) < 2:
            return sequences
        elif len(points) == 2:
            return [[inj] for inj in injections(points[0],points[1])]
        else:
            for inj in injections(points[0],points[1]):
                coker = set(points[1]).difference(inj.values())
                sequences.extend([[inj]+sequence for sequence in CTMinus.enumerate_gens_helper([list(coker)]+points[2:])])
        return sequences
