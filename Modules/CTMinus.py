from __future__ import annotations
from Tangles.Tangle import *
from Modules.Bimodule import *
from Tangles.Functions import *
import copy


# represents a pair of partial bijection overlaid
class StrandDiagram:
    def __init__(self, etangle: ETangle, left_strands: Dict, right_strands: Dict):
        self.etangle = etangle
        self.left_strands = left_strands
        self.right_strands = right_strands


def type_da(etangle: ETangle) -> Bimodule:
    left_algebra = AMinus(etangle.left_signs())
    right_algebra = AMinus(etangle.right_signs())
    gens = enumerate_gens([etangle.left_points(), etangle.middle_points(), etangle.right_points()])
    maps = None  # TODO

    return Bimodule(left_algebra, right_algebra, gens, maps)

def dplus(sd: StrandDiagram):
    strands = sd.right_strands
    keys = strands.keys()
    out = {}
    for key1 in keys:
        for key2 in keys:
            if key2 < key1 and strands[key2]>strands[key1]:
                res = resolveplus(sd,key1,key2)
                if res[1] != 0 and res[0] in out.keys():
                    out[res[0]] = res[1]+out[res[0]]
                elif res[1] != 0:
                    out[res[0]] = res[1]
    return out

def resolveplus(sd: StrandDiagram, i, j):
    strands = sd.right_strands
    # if double crossing black, return none
    for s in strands.keys() & set(range(j,i)):
        if strands[i]<strands[s]<strands[j]: 
            return [None,0]
        
    # output
    out = StrandDiagram(sd.etangle,copy.deepcopy(sd.left_strands),copy.deepcopy(sd.right_strands))
    out.right_strands[j] = strands[i]
    out.right_strands[i] = strands[j]


    # calculate coefficient from orange correctly
    pos = sd.etangle.position
    if sd.etangle.etype == ETangle.Type.CAP:
        checkrange = range(max(strands[i],j),min(i,strands[j]))
        c = sd.etangle.polyring.one()
        for k in checkrange:
            index = k
            if k >= pos - 1: index = k+2
            if sd.etangle.signs[index] == -1: 
                return [None,0]
            else:
                c = c*sd.etangle.polyring['U'+str(sd.etangle.middle[index]+1)]
        return [out,c]
    elif sd.etangle.etype == ETangle.Type.CUP:
        checkrange = range(max(strands[i],j+j>=pos),min(i+i>=pos,strands[j]))
        c = sd.etangle.polyring.one()
        for k in checkrange:
            if sd.etangle.signs[k] == -1:
                return [None,0]
            else:
                c = c*sd.etangle.polyring['U'+str(sd.etangle.middle[k]+1)]
        return [out,c]
    else:
        checkrange = range(max(strands[i],j),min(i,strands[j]))
        c = sd.etangle.polyring.one()
        for k in checkrange:
            if sd.etangle.signs[k] == -1:
                return [None,0]
            else:
                c = c*sd.etangle.polyring['U'+str(sd.etangle.middle[k]+1)]
        return [out,c]


# points - a list of sets of points
def enumerate_gens(points):
    def reverse_injection(d):
        return {v: k for k, v in d.items()}

    sequences = []
    if len(points) < 2:
        return sequences
    elif len(points) == 2:
        return [[pb] for pb in partial_bijections(points[0], points[1])]
    else:
        mid = len(points) // 2
        for pb in partial_bijections(points[mid], points[mid + 1]):
            r, s = list(reversed(points[:mid + 1])), points[mid + 1:]

            r[0] = list(set(r[0]).difference(pb.keys()))
            s[0] = list(set(s[0]).difference(pb.values()))

            t = list(enumerate_gens_helper(s))
            sequences.extend([[reverse_injection(inj) for inj in reversed(s1)]
                              + [pb] + s2 for s1 in enumerate_gens_helper(r) for s2 in t])
        return sequences


def enumerate_gens_helper(points):
    if len(points) == 2:
        for inj in injections(points[0], points[1]):
            yield [inj]
    elif len(points) > 2:
        for inj in injections(points[0], points[1]):
            coker = set(points[1]).difference(inj.values())
            for sequence in enumerate_gens_helper([list(coker)] + points[2:]):
                yield [inj] + sequence
