from __future__ import annotations
from typing import List
from Tangles.Tangle import *
from Modules.Bimodule import *
from Tangles.Functions import *
from Modules.StrandDiagram import *
import copy


def type_da(etangle: ETangle) -> TypeDA:
    strand_diagrams = [StrandDiagram(etangle, left_strands, right_strands)
                       for left_strands, right_strands in
                       enumerate_gens([etangle.left_points(), etangle.middle_points(), etangle.right_points()])]
    maps = sum((delta1_1(x) for x in strand_diagrams), []) + \
        [delta1_2(x, a) for x in strand_diagrams
         for a in etangle.right_algebra.left_gens(list(x.left_strands.keys()))]

    return TypeDA.from_strand_diagrams(etangle.left_algebra, etangle.right_algebra, strand_diagrams, maps)


def delta1_1(x: StrandDiagram) -> List[Bimodule.Edge]:
    out = []
    out += [Bimodule.Edge(x, y, c, (x.left_idempotent(),), tuple()) for y, c in dplus(x).items()]
    out += [Bimodule.Edge(x, y, c, (x.left_idempotent(),), tuple()) for y, c in dminus(x).items()]
    # out += [Bimodule.Edge(x, y, c, (x.left_idempotent(),), tuple()) for y, c in dmixed(x).items()]
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
    t = sd.etangle.type

    # if double crossing black, return none
    for s in strands.keys() & set(range(j,i)):
        if strands[i]<strands[s]<strands[j]:
            return [None,zero]

    # output
    left_strands = dict(sd.left_strands)
    right_strands = dict(sd.right_strands)
    right_strands[j] = strands[i]
    right_strands[i] = strands[j]
    out = StrandDiagram(sd.etangle, left_strands, right_strands)

    # calculate coefficient from orange correctly
    pos = sd.etangle.position
    c = sd.etangle.polyring.one()

    if t == ETangle.Type.UNDER:
        checkrange = range(max(strands[i],j),min(i,strands[j]))
        signs = sd.etangle.right_signs()
        for k in checkrange:
            if signs[k] == -1:
                return [None,zero]
            else:
                c = c*sd.etangle.polyring['U'+str(k+1)]
        return [out,c]
    elif t == ETangle.Type.OVER:
        checkrange = range(max(strands[i],j+j>=pos),min(i,strands[j]+strands[j]>=pos))
        signs = sd.etangle.left_signs()
        for k in checkrange:
            if signs[k] == -1:
                return [None,0]
            else:
                c = c*sd.etangle.polyring['U'+str(k+1)]
        return [out,c]
    elif t == ETangle.Type.CAP:
        checkrange = range(max(j - j>=pos,strands[i]), min(i-i>=pos,strands[j]))
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
    elif t == ETangle.Type.CUP:
        signs = sd.etangle.right_signs()
        checkrange = range(max(strands[i],j+j>=pos),min(i+i>=pos,strands[j]))
        for k in checkrange: 
            if sd.etangle.signs[k] == -1:
                return [None,zero]
            else: 
                c = c*sd.etangle.polyring['U'+str(k+1)]
        return [out,c]


def resolveminus(sd:StrandDiagram, i, j):
    zero = sd.etangle.polyring.zero()
    strands = sd.left_strands
    signs = sd.etangle.left_signs()

    #check if we need to double cross black to introduce crossing
    for s in strands.keys() & set(range(j,i)):
        if strands[j]<strands[s]<strands[i]:
            return [None,zero]

    #output
    left_strands = dict(sd.left_strands)
    right_strands = dict(sd.right_strands)
    left_strands[j] = strands[i]
    left_strands[i] = strands[j]
    out = StrandDiagram(sd.etangle, left_strands, right_strands)

    # calculate coefficient from orange
    pos = sd.etangle.position
    c = sd.etangle.polyring.one()
    t = sd.etangle.etype

    if t == ETangle.Type.OVER:
        checkrange = range(max(j,strands[j]),min(i,strands[i]))
        for k in checkrange:
            if signs[k] == 1:
                return [None,zero]
            else:
                c = c*sd.etangle.polyring['U'+str(k+1)]
        return [out,c]
    elif t == ETangle.Type.UNDER:
        checkrange = range(max(j+(j>=pos),strands[j]),min(i,strands[i]+(strands[i]>=pos)))
        for k in checkrange:
            if k not in [pos-1,pos]:
                if signs[k] == 1:
                    return [None,zero]
                else:
                    c = c*sd.etangle.polyring['U'+str(k+1)]
            if k == pos:
                if strands[j] < k:
                    if signs[k] == 1:
                        return [None,zero]
                    else: 
                        c = c*sd.etangle.polyring['U' + str(k-1+1)]
            if k == pos-1:
                if strands[i] > pos:
                    if signs[k] == 1:
                        return [None,zero]
                    else:
                        c = c*sd.etangle.polyring['U'+str(k+1+1)]
        return [out,c]
    elif t == ETangle.Type.CAP:
        checkrange = range(max(j+(j>=pos),strands[j]),min(i+(i>=pos),strands[i]))
        for k in checkrange:
            if signs[k] == 1:
                return [None,zero]
            else:
                c = c*sd.etangle.polyring['U'+str(k+1)]
        return [out,c]
    elif t == ETangle.Type.CUP:
        checkrange = range(max(j - (j>=pos),strands[j]),min(i-(i>=pos),strands[i]))
        signs = sd.etangle.right_signs()
        for k in checkrange:
            if k>=pos - 1:
                if signs[k+2] == 1:
                    return [None,zero]
                else: 
                    c = c*sd.etangle.polyring['U'+str(k+2+1)]
            else:
                if signs[k] == 1:
                    return [None,zero]
                else:
                    c = c*sd.etangle.polyring['U'+str(k+1)]
        return [out,c]


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
