

class StrandDiagram:
    def __init__(self, parent, bijection):
        assert (parent.degree + 1) == len(bijection)
        # need to assert it is a valid partial bijection:
        # 1) no number is repeated twice
        # 2) every number is positive and < len(bijection)

        # parent is a sign sequence
        self.parent = parent
        self.degree = parent.degree
        # bijection is a list of size degree + 1. bijection[i] says where i maps to, 
        # or -1 if there is no strand out of i
        self.bijection = bijection

    def parent(self):
        return self.parent

    def degree(self):
        return self.degree
    
    def bijection(self):
        return self.bijection

    def __str__(self):
        return str(self.bijection)+", "+str(self.parent)

    def __repr__(self):
        return str(self.bijection)+", "+str(self.parent)

    # THIS IS NOT IMPLEMENTED
    def __mul__(self, other):
        for i in self.bijection: 
            if (i != -1) and (other.bijection[i] == -1): return 0

    def __eq__(self, other):
        if (self.parent == other.parent) and (self.bijection == other.bijection): return true
        else: return false

    # THIS IS BAD BAD BAD
    def __hash__(self):
        return hash(str(self))


    def differential(self):
        bijection = self.bijection
        diff = SumOfStrand({},self.parent)
        for i in range(1,len(bijection)):
            if bijection[i] != -1:
                for j in range(0,i):
                    if bijection[j]>bijection[i]:
                        resolution = resolve(self,i,j)
                        diff.add(resolution)
        return diff

def resolve(sd,i,j):
    bij = sd.bijection
    parent = sd.parent

    
    checkrange = range(max(bij[i],j),min(i,bij[j]))
    
    # check if black strands double cross
    for k in range(j,i): 
        if (bij[i]<bij[k]) and (bij[k]<bij[j]):
            return [sd,0]

    newgen = list(bij)
    newgen[i]=bij[j]
    newgen[j]=bij[i]
    newgen = StrandDiagram(parent,newgen)

    # if it crosses no orange strands return coeff 1
    if checkrange == []:
        return [newgen,1]
    
    # if a black strand pulls across a negative (left) oriented orange strand, don't count it
    for k in checkrange:
        if parent.sequence[k] == -1:
            return [sd,0]
    return [newgen,checkrange]

