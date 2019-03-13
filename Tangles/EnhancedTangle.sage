#load("AlgElement.sage")

import copy

#An EnhancedETangle is an ETangle equipped with a (left) partial bijection bijminus and a (right) partial bijplus
#where bijminus is a partial bijection from the set of points on the left to the set of points in the middle,
#and bijplus is a partial bijection from the set of points in the middle to the set of points on the right
#such that the set of middle points is the disjoint union of the image of bijminus and the domain of bijplus
#Basically, we are overlaying black strands on orange strands
#Should be built with ETangle as parent, but we haven't configured it that way

class EnhancedETangle():

    def __init__(self, parent, bijminus, bijplus):
        # Parent is an ETangle,
        #bijminus is a dictionary taking in left points and returning middle points
        #bijplus is a dictionary taking in middle points and returning right points
        # Indexing on bijections starts on 0 as the bottom strand

        #These are the images and domains of our partial bijections--note that they are unordered
        self.minusimage = set(bijminus[a] for a in bijminus)
        self.plusimage = set(bijplus[a] for a in bijplus)
        self.minusdomain = set(bijminus)
        self.plusdomain = set(bijplus)

        #We check that the disjoint union of minusimage and plusdomain is the middle points
        assert self.minusimage & self.plusdomain == set([])
        assert self.minusimage | self.plusdomain == set(parent.middle_points())

        # need to assert this is a valid bijection still
        # i.e. max value is at most parent.degree + 1, and
        # no value is hit twice
        self.parent = parent
        self.position = parent.position
        self.signs = parent.signs

        self.bijplus = bijplus
        self.bijminus = bijminus
        self.degree = parent.degree
        self.etype = parent.etype

    # Access Methods

    #The ETangle of our object
    def parent(self):
        return self.parent

    def etangle(self):
       return self.parent

    #The sign sequence of the ETangle of our object
    def signs(self):
        return self.signs

    #The position where the orange strands of ETangle are doing stuff
    def position(self):
      return self.position

    #The degree of the sign sequence of the ETangle of our object
    def degree(self):
        return self.degree

    def len(self):
        return self.degree

    def bijminus(self):
        return self.bijminus

    def bijplus(self):
        return self.bijplus

    #Takes in an enhanced elementary tangle and produces a dictionary with keys that are enhanced elementary Tangles
    #and outputs their coefficients in F2[U1 ... U_npos]
    #Note that this is a formal sum of elements in the algebra, and no longer an (enhanced) elementary tangle
    def Dplus(self):
        npos = self.degree
        var_names = ['U%s'%p for p in range(1,npos + 1)]
        #Polynomial ring where all of our algebra takes place
        #Ideally, this ring would exist globally and we wouldn't need to coerce things in this function, but I don't know how to make that happen
        R = PolynomialRing(GF(2),npos,var_names)
        #Output is the image of Dplus. We will adjoin elements to it later
        output = {}
        #the right sign sequence--that is, the sign sequence on the right side of our elementary tangle
        #Also, it is the right sign sequence for us to be using
        rightsigns = list(self.signs)
        #Because signs refers to different things for different etypes, we define rightsigns in cases
        if self.etype == ETangle.Type.CAP:
            rightsigns[self.position - 1] = 0
            rightsigns[self.position] = 0
        if self.etype in (ETangle.Type.OVER, ETangle.Type.UNDER):
            signpos = rightsigns[self.position]
            signposneg = rightsigns[self.position - 1]
            rightsigns[self.position - 1] = signpos
            rightsigns[self.position] = signposneg
        #We iterate over all i, j in the domain of bijplus, looking for black-black crossings to uncross
        for i in self.plusdomain:
            for j in self.plusdomain:
                #If we see a black-black crossing, we uncross it
                #We name the new bijection newgen and the associated elementary tangle newtangle
                #newtangle is a candidate term in our formal sum; we'll determine the coefficients later
                if i > j and self.bijplus[i] < self.bijplus[j]:
                    #We want a deep copy so we don't modify bijplus in place
                    newgen = copy.deepcopy(self.bijplus)
                    newgen[i]=self.bijplus[j]
                    newgen[j]=self.bijplus[i]
                    #The tangle we will be looking for the coefficient of
                    newtangle = EnhancedETangle(self.parent,self.bijminus,newgen)
                    #check if black strands double-cross in newtangle: if they do, we don't include this term in the output
                    doublecross = False
                    for k in xrange(j+1,i):
                        try:
                            if (self.bijplus[i] < self.bijplus[k] and newgen[i] > newgen[k]) or (self.bijplus[j] > self.bijplus[k] and newgen[j] < newgen[k]):
                                doublecross = True
                                break
                        except:
                            pass
                    if doublecross == False:
                        #ODC = OrangeDoubleCrossings. Keeps track of where (positive or negative) orange strands are double crossed
                        #The 0th entry of ODC will be the coefficient of the term. If it is 0, there is no need to include the term
                        ODC = [1] + [0]*npos
                        #The only place where new double crossings could occur is in the range that was changed, namely from j to i
                        for k in xrange(j+1,i):
                            #Orange double crossings are dealt with in cases
                            if self.etype == ETangle.Type.OVER:
                                #OVER has complicated double-crossing relations around the position
                                if k in (self.position,self.position + 1):
                                    if (j < self.position and self.bijplus[j] > self.position and
                                    i > self.position and self.bijplus[i] <= self.position):
                                        ODC[self.position] = 1
                                    elif (i > self.position and self.bijplus[i] < self.position and
                                    j < self.position and self.bijplus[j] >= self.position):
                                        ODC[self.position+1] = 1
                                else:
                                    if (self.bijplus[i] < k and newgen[i] >= k) or (self.bijplus[j] >= k and newgen[j] < k):
                                        ODC[k] = 1
                            elif self.etype == ETangle.Type.UNDER:
                                #UNDER is less complicated but we need to keep track of which Uk was crossed
                                if (self.bijplus[i] < k and newgen[i] >= k) or (self.bijplus[j] >= k and newgen[j] < k):
                                    if k == self.position:
                                        ODC[self.position+1] = 1
                                    elif k == self.position + 1:
                                        ODC[self.position] = 1
                                    else:
                                        ODC[k] = 1
                            else:
                                #In any other case, checking for double-crossings is straightforward
                                if (self.bijplus[i] < k and newgen[i] >= k) or (self.bijplus[j] >= k and newgen[j] < k):
                                    ODC[k] = 1
                        #We now generate our coefficients
                        for ell in xrange(1,len(ODC)):
                            #If any of the double crossings are negatively oriented, we don't need this term in our sum
                            if ODC[ell] == 1 and rightsigns[ell-1] == -1:
                                ODC[0] = 0
                                break
                            #For each positively oriented double crossed orange strand, we will have a Uk in our coefficient
                            elif ODC[ell] == 1 and rightsigns[ell-1] == 1:
                                if self.etype == ETangle.Type.CUP and k in (self.position, self.position + 1):
                                    a = self.position
                                    #This multiplies ODC[0] by U_{self.position}
                                    ODC[0] *= R('U%s'%a)
                                else:
                                    #This multiplies ODC[0] by U_k
                                    ODC[0] *= R('U%s'%ell) #We coerce Uk to live in the ring GF(2)[U1 ... Uk]
                        if ODC[0] != 0:
                            output[newtangle] = ODC[0]
                        print ODC
        return output
