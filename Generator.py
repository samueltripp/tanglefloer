

class Generator:
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


    # ACCESS METHODS
    def parent(self):
        return self.parent

    def degree(self):
        return self.degree
    
    def bijection(self):
        return self.bijection

    # UNDERLYING METHODS
    def __str__(self):
        return str(self.bijection)+", "+str(self.parent)

    def __repr__(self):
        return str(self.bijection)+", "+str(self.parent)

    def __eq__(self, other):
        if (self.parent == other.parent) and (self.bijection == other.bijection): return true
        else: return false

    # I think this is bad practice but don't know. It does 
    # uniquely identify elements, because two Generators are equal
    # if and only if their string reps are
    def __hash__(self):
        return hash(str(self))


    # ARITHMETIC METHODS
    # This adds two elements
    # Takes in two Generators and returns their product, an AlgElement
    def __add__(self,other):
        assert self.parent == other.parent
        sum = AlgElement({},self.parent)
        sum.add([self,1])
        sum.add([other,1])
        return sum

    # This multiplies two elements
    # Takes in two Generators and returns their product, an AlgElement
    def __mul__(self, other):
        # If these generators aren't even for the same SignSequence, throw an error
        assert self.parent == other.parent
        product = AlgElement({},self.parent)

        # Return 0 AlgElement if the ends of self strands don't match the beginnings of other strands
        if (not ends_match(self,other)): return product

        bijone = self.bijection
        bijtwo = other.bijection

        # If two black strings double cross, return 0 AlgElement
        for i in range(1,len(bijone)):
            if bijone[i] != -1:
                for j in range(0,i):
                    if (bijone[j] > bijone[i]) and (bijtwo[bijone[j]] < bijtwo[bijone[i]]): return product

        # coeffdev is a list that counts how many times each orange strand is double crosswed
        coeffdev = [0]*(len(bijone)-1)
        for i in range(0,len(bijone)):
            if bijone[i] != -1:
                if bijone[i] > i:
                    checkrange = range(max(i,bijtwo[bijone[i]]),bijone[i])
                else: checkrange = range(bijone[i],min(i,bijtwo[bijone[i]]))
                for k in checkrange:
                    coeffdev[k] = coeffdev[k] + 1

        # This is creating the new generator that is just the composition of self and other
        newbij = [-1]*len(bijone)
        for i in range(0,len(bijone)):
            if bijone[i] != -1:
                newbij[i] = bijtwo[bijone[i]]
        newgen = Generator(self.parent,newbij)

        # Add this new generator with appropriate coefficients to the AlgElement, and return
        product.add([newgen,coeffdev])
        return product



    # Returns an AlgElement of the differential of this Generator
    # Returns the zero element if your generator has no crossings
    def differential(self):
        bijection = self.bijection
        diff = AlgElement({},self.parent)
        for i in range(1,len(bijection)):
            if bijection[i] != -1:
                for j in range(0,i):
                    if bijection[j]>bijection[i]:
                        resolution = resolve(self,i,j)
                        diff.add(resolution)
        return diff

    # Use this method to add or multiply Generators to AlgElements
    def toAlgElement(self):
        return AlgElement({self:1},self.parent)

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
    newgen = Generator(parent,newgen)

    # if it crosses no orange strands return coeff 1
    if checkrange == []:
        return [newgen,1]
    
    # if a black strand pulls across a negative (left) oriented orange strand, don't count it
    for k in checkrange:
        if parent.sequence[k] == -1:
            return [sd,0]

    returnrange = [0]*len(sd.parent.sequence)
    for k in checkrange:
        returnrange[k] = 1
    return [newgen,returnrange]



def ends_match(sd1,sd2):
    bijlength = len(sd1.bijection)

    selfright = [-1]*bijlength
    otherleft = list(selfright)

    for i in range(0,bijlength):
        if (sd1.bijection[i] != -1):
            selfright[sd1.bijection[i]] = 1
            
    for i in range(0,bijlength):
        if (sd2.bijection[i] != -1):
            otherleft[i] = 1        
    return selfright == otherleft

    
