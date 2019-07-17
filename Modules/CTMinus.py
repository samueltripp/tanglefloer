from __future__ import annotations
from typing import List
from Tangles.Tangle import *
from Modules.Bimodule import *
from Tangles.Functions import *
from Modules.StrandDiagrams import *
import copy


def type_da(etangle: ETangle) -> Bimodule:
    strand_diagrams = [StrandDiagram(etangle, left_strands, right_strands)
                       for left_strands, right_strands in
                       enumerate_gens([etangle.left_points(), etangle.middle_points(), etangle.right_points()])]
    maps = sum((delta1_1(x) for x in strand_diagrams), []) + \
        [delta1_2(x, a) for x in strand_diagrams
         for a in etangle.right_algebra.left_gens(list(x.left_strands.keys()))]

    return Bimodule.from_strand_diagrams(etangle.left_algebra, etangle.right_algebra, strand_diagrams, maps)


def delta1_1(x: StrandDiagram) -> List[Bimodule.Edge]:
    out = []
    out += [Bimodule.Edge(x, e.target_diagram, e.c, (x.left_idempotent(),), tuple()) for e in dplus(x)]
    # out += [Bimodule.Edge(x, e.target_diagram, e.c, (x.left_idempotent(),), tuple()) for e in dminus(x)]
    # out += [Bimodule.Edge(x, e.target_diagram, e.c, (x.left_idempotent(),), tuple()) for e in dmixed(x)]
    out += [deltal(x)]
    return out


def delta1_2(x: StrandDiagram, a: AMinusElement) -> Bimodule.Edge:
    e = m2(x, a)
    return Bimodule.Edge(x, e.target_diagram, e.c, (x.left_idempotent(),), (a,))


def m2(x: StrandDiagram, a: AMinusElement) -> Bimodule.Edge:
    pass  # TODO


def deltal(x: StrandDiagram) -> Bimodule.Edge:
    pass  # TODO


def dplus(sd: StrandDiagram):
    strands = sd.right_strands
    keys = strands.keys()
    zero = sd.etangle.polyring.zero()
    out = {}
    for key1 in keys:
        for key2 in keys:
            if key2 < key1 and strands[key2]>strands[key1]:
                res = resolveplus(sd,key1,key2)
                if res[1] != zero and res[0] in out.keys():
                    out[res[0]] = res[1]+out[res[0]]
                elif res[1] != zero:
                    out[res[0]] = res[1]
    return out

def dminus(sd: StrandDiagram):
    strands = sd.left_strands
    keys = strands.keys()
    out = {}
    zero = sd.etangle.polyring.zero()
    for key1 in keys:
        for key2 in keys:
            if key2 < key1 and strands[key2]<strands[key1]:
                res = resolveminus(sd,key1,key2)
                if res[1] != zero and res[0] in out.keys():
                    out[res[0]] = res[1]+out[res[0]]
                elif res[1] != zero:
                    out[res[0]] = res[1]
    return out

def resolveplus(sd: StrandDiagram, i, j):
    zero = sd.etangle.polyring.zero()
    strands = sd.right_strands
    # if double crossing black, return none
    for s in strands.keys() & set(range(j,i)):
        if strands[i]<strands[s]<strands[j]:
            return [None,zero]

    # output
    out = StrandDiagram(sd.etangle,copy.deepcopy(sd.left_strands),copy.deepcopy(sd.right_strands))
    out.right_strands[j] = strands[i]
    out.right_strands[i] = strands[j]

    # calculate coefficient from orange correctly
    pos = sd.etangle.position
    if sd.etangle.etype == ETangle.Type.UNDER:
        checkrange = range(max(strands[i],j),min(i,strands[j]))
        c = sd.etangle.polyring.one()
        signs = sd.etangle.right_signs()
        for k in checkrange:
            if signs[k] == -1:
                return [None,zero]
            else:
                c = c*sd.etangle.polyring['U'+str(k+1)]
        return [out,c]
    elif sd.etangle.etype == ETangle.Type.OVER:
        checkrange = range(max(strands[i],j+j>=pos),min(i,strands[j]+strands[j]>=pos))
        c = sd.etangle.polyring.one()
        signs = sd.etangle.left_signs()
        for k in checkrange:
            if signs[k] == -1:
                return [None,0]
            else:
                c = c*sd.etangle.polyring['U'+str(k+1)]
        return [out,c]
    elif sd.etangle.etype == ETangle.Type.CAP:
        checkrange = range(max(j - j>=pos,strands[i]), min(i-i>=pos,strands[j]))
        c = sd.etangle.polyring.one()
        signs = sd.etangle.left_signs()
        for k in checkrange:
            if k>=pos - 1:
                if signs[k+2] == -1:
                    return [None,zero]
                else: c=c*sd.etangle.polyring['U'+str(k+2+1)]
            else: 
                if signs[k] == -1:
                    return [None,zero]
                else: c = c*sd.etangle.polyring['U'+str(k+1)]
        return [out,c]
    elif sd.etangle.etype == ETangle.Type.CUP:
        checkrange = range(max(strands[i],j+j>=pos),min(i+i>=pos,strands[j]))
        c = sd.etangle.polyring.one()
        for k in checkrange: 
            if sd.etangle.signs[k] == -1:
                return [None,zero]
            else: 
                c = c*sd.etangle.polyring['U'+str(k+1)]
        return [out,c]


def resolveminus(sd:StrandDiagram, i, j):
    zero = sd.etangle.polyring.zero()
    strands = sd.right_strands

    #check if we need to double cross black to introduce crossing
    
    



# points - a list of sets of points
def enumerate_gens(points):
    sequences = []
    if len(points) < 2:
        return sequences
    elif len(points) == 2:
        return [[pb] for pb in partial_bijections(points[0], points[1])]
    else:
        for pb in partial_bijections(points[0], points[1]):
            coker = set(points[1]).difference(pb.values())
            sequences.extend(
                [[pb] + sequence for sequence in enumerate_gens_helper([list(coker)] + points[2:])])
    return sequences


def enumerate_gens_helper(points):
    sequences = []
    if len(points) < 2:
        return sequences
    elif len(points) == 2:
        return [[inj] for inj in injections(points[0], points[1])]
    else:
        for inj in injections(points[0], points[1]):
            coker = set(points[1]).difference(inj.values())
            sequences.extend(
                [[inj] + sequence for sequence in enumerate_gens_helper([list(coker)] + points[2:])])
    return sequences
