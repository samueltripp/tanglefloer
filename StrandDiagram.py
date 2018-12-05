

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

def differential(sd):
    bijection = sd.bijection
    diff = SumOfStrand({})
    for i in range(1,len(bijection)):
        if bijection[i] != -1:
            for j in range(0,i):
                if bijection[j]>bijection[i]:
                    diff.add(resolve(sd,i,j),1)
    return diff

def resolve(sd,i,j):
    range = range(max(sd.bijection[i],j),min(i,sd.bijection[j]))

    # if a black strand pulls across a negative (left) oriented orange strand, don't count it
    for k in range:
        if sd.parent.sequence[k] == -1:
            return 0
        
    # if a black strand double crosses another, don't count it
    for k in range(j,i):
        if (sd.bijection[i]<sd.bijection[k]) and (sd.bijection[k]<sd.bijection[j]):
            return 0

    # if neither of these happen, find the right coefficients and return them
    return sd
