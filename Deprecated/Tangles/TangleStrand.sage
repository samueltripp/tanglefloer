class ETangleStrand():

    def __init__(self, etangle, bijminus,bijplus):
        self.etangle = etangle
        self.bijminus = bijminus
        self.bijplus = bijplus

        self.minusimage = set(bijminus[i] for i in bijminus)
        self.minusdomain = set(bijminus)
        self.plusimage = set(bijplus[i] for i in bijplus)
        self.plusdomain = set(bijplus)

        # Assert a valid pair of bijections
        assert self.validPair()
        
        self.degree = etangle.degree
        self.etype = etangle.etype
        self.position = etangle.position
        self.signs = etangle.signs
        self.ring = etangle.ring
        


    def dplus(self):
        # A dictionary between etanglestrands and coefficients for output
        output = {}

        rightsigns = self.etangle.right_signs()
        for i in self.plusdomain:
            for j in self.plusdomain:
                if i > j and bijplus[i] < bijplus[j]:
                    resolution = dplusresolve(self,i,j)
                    if resolution[1] != 0:
                        if resolution[0] in output: 
                            output[resolution[0]] = output[resolution[0]]+resolution[1]
                        else output[resolution[0]] = resolution[1]

    def dplusresolve(i,j):
        


    def validPair(self):
        return self.minusimage & self.plusdomain == set([]) and self.minusimage | self.plusdomain == set(self.etangle.middle_points())
