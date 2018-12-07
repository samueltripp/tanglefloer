class AlgElement:
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
                coeff = self.find_coeff(range)
            if gen in dict:
                dict[gen] = dict[gen]+coeff
            else: 
                dict[gen] = coeff

    def find_coeff(self,inrange):
        gens = self.gens
        poscount = self.poscount
        retcoeff = 1
        if all([v == 0 for v in inrange]): return 1
        for i in range(0,len(inrange)):
            if (inrange[i]>0) and (poscount[i] == -1):
                return 0
            else: 
                retcoeff = retcoeff * gens[poscount[i]-1]**inrange[i]

        return retcoeff
