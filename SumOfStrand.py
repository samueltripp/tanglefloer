class SumOfStrand:
    def __init__(self,dict,parent):
        self.dict = dict
        self.parent = parent
        self.poscount = parent.poscount
        self.polyring = self.initpolyring(parent)
        self.gens = self.polyring.gens()

    def dict(self):
        return self.dict

    def parent(self):
        return self.parent

    def polyring(self):
        return self.polyring

    def __eq__(self,other):
        if (self.dict == other.dict) and (self.parent == other.parent): return true
        else: return false

    @staticmethod
    def initpolyring(parent):
        npos = parent.npos
        var_names = ['U%s'%p for p in range(1,npos + 1)]
        return PolynomialRing(GF(2),npos,var_names)

    def add(self,input):
        dict = self.dict
        gen = input[0]
        range = input[1]
        gens = self.gens
        poscount = self.poscount

        if range != 0:
            coeff = 1
            if range != 1:
                for k in range:
                    coeff = coeff*gens[poscount[k]-1]
            if gen in dict:
                dict[gen] = dict[gen]+coeff
            else: 
                dict[gen] = coeff


    def pop(self,key):
        self.dict.pop(key,None)
