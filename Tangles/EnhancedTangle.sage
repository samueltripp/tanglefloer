#An EnhancedETangle is an ETangle equipped with a (left) partial bijection bijminus and a (right) partial bijplus
#where bijminus is a partial bijection from the set of points on the left to the set of points in the middle,
#and bijplus is a partial bijection from the set of points in the middle to the set of points on the right
#such that the set of middle points is the disjoint union of the image of bijminus and the domain of bijplus
#Basically, we are overlaying black strands on orange strands
class EnhancedETangle(ETangle):

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
        ETangle.__init__(self, parent.etype, parent.signs, parent.position)
        self.parent = parent
        self.bijplus = bijplus
        self.bijminus = bijminus

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

    #This code is not actually Dplus yet
    #Instead, it takes in self and returns a list of copies of self, each of which has one crossing uncrossed
    def Dplus(self):
        output = {}
        for i, j in self.plusdomain:
            #wherever there are crossings, we uncross and create a new term
            if i > j and bijplus[i] < bijplus[j]:
                newgen = copy.deepcopy(bijplus)
                newgen[i]=bijplus[j]
                newgen[j]=bijplus[i]
                #check if black strands double-cross: if they do, we don't include this term in the output
                doublecross = False
                for k in range(j,i):
                    try:
                        if (bij[i]<bij[k]) and (bij[k]<bij[j]):
                            doublecross = True
                            break
                    except:
                        pass
                if doublecross == False:
                    #if we doublecross left-oriented orange strands, we don't include this term in the output
                    if i < bijplus[i] and newgen[j] < bijplus[j]:
                        orangestrandcrossing == False
                        for k in [i-1,newgen[i]]:
                            if self.signs[k] == -1:
                                orangestrandcrossing == True
                                break
                        #SOMETHING GOES HERE
                    if i > bijplus[i] and newgen[j] > bijplus[j]:
                        orangestrandcrossing = False
                        for k in [i-1,newgen[i]]:
                            if self.signs[k] == -1:
                                orangestrandcrossing == True
                                break
                        #SOMETHING GOES HERE
        return output
