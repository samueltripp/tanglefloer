

class StrandDiagram:
    def __init__(self, parent, bijection):
        assert (parent.degree + 1) == len(bijection)
        # need to assert it is a valid partial bijection:
        # 1) no number is repeated twice
        # 2) every number is positive and < len(bijection)
        self.parent = parent
        self.degree = parent.degree
        self.bijection = bijection

    def parent(self):
        return self.parent

    def degree(self):
        return self.degree
    
    def bijection(self):
        return self.bijection

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
    return sd
                    
