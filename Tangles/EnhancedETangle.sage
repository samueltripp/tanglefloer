#An EnhancedETangle is an ETangle equipped with a (left) partial bijection bijminus and a (right) partial bijplus
#where bijminus is a partial bijection from the set of points on the left to the set of points in the middle,
#and bijplus is a partial bijection from the set of points in the middle to the set of points on the right
#such that the set of middle points is the disjoint union of the image of bijminus and the domain of bijplus
#Basically, we are overlaying black strands on orange strands
class EnhancedETangle:

    def __init__(self, parent, bijminus, bijplus):
        # Parent is an ETangle,
        #bijminus is a dictionary taking in left points and returning middle points
        #bijplus is a dictionary taking in middle points and returning right points
        # Indexing on bijections starts on 0 as the bottom strand

        #These are the images and domains of our partial bijections--note that they are unordered
        minusimage = set(bijminus[a] for a in bijminus)
        plusimage = set(bijplus[a] for a in bijplus)
        minusdomain = set(bijminus)
        plusdomain = set(bijplus)

        #We check that the disjoint union of minusimage and plusdomain is the middle points
        assert minusimage & plusdomain == set([])
        assert minusimage | plusdomain == set(parent.middle_points())

        # need to assert this is a valid bijection still
        # i.e. max value is at most parent.degree + 1, and
        # no value is hit twice
        self.parent = parent
        self.signs = parent.signs
        self.degree = len(parent.signs)
        self.position = parent.position

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

    def Dplus(self):
      
