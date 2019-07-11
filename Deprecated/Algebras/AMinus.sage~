# class to represent A^-(P) for some sign sequence P
class AMinus:
    def __init__(self, sign_sequence):
        self.sign_sequence = sign_sequence
        self.signSequence = SignSequence(list(sign_sequence))
        self.positives = tuple(i for i in range(len(sign_sequence)) if sign_sequence[i] == 1)
        self.gens = self.enumerate_gens(len(self.sign_sequence)+1)
        self.polyring = self.initpolyring()
        self.algebra = self.initalgebra()

    def initpolyring(self):
        if len(self.positives) == 0:
            return GF(2)
        else: 
            return PolynomialRing(GF(2),len(self.positives),['U%s'%p for p in range(1,len(self.positives)+1)])

    def initalgebra(self):
        mats = [self.matrix_representation(i) for i in range(len(self.gens))]
        return FiniteDimensionalAlgebra(self.polyring,mats)

        
    # populates self.gens with every generator
    def enumerate_gens(self, a):
        gentuples =  self.enumerate_gens_helper(a, list(range(a)))
        gens = [Generator(self.signSequence,list(gentuples[i])) for i in range(len(gentuples))]
        return gens

    # helps
    def enumerate_gens_helper(self, a, codomain):
        if a == 0:
            return [()]
        if len(codomain) == 0:
            return [(-1,)*a]
        gens = [l+(-1,) for l in self.enumerate_gens_helper(a-1,codomain)]
        for i in range(len(codomain)):
            gens.extend([l+(codomain[i],) for l in self.enumerate_gens_helper(a-1,codomain[:i]+codomain[i+1:])])
        return gens

    # multiplies two generators
    def multiply(self, gena, genb):
        gen1 = gena.bijection
        gen2 = genb.bijection
        # if ends of strands don't match
        for i in range(len(gen1)):
            if (i not in gen1.values()) != (gen2[i] == -1):
                return (0,None)

        # if two black strands double cross
        for i in range(1,len(gen1)):
            if gen1[i] != -1:
                for j in range(i):
                    if (gen1[j] > gen1[i]) and (gen2[gen1[j]] < gen2[gen1[i]]):
                        return (0,None)

        # counts how many times each orange strand is double-crossed
        coefficient = 1
        for i in range(0,len(gen1)):
                if gen1[i] != -1:
                    if gen1[i] > i:
                        checkrange = range(max(i,gen2[gen1[i]]),gen1[i])
                    else: checkrange = range(gen1[i],min(i,gen2[gen1[i]]))
                    for k in checkrange:
                        if k not in self.positives:
                            return (0,None)
                        else:
                            coefficient = coefficient*self.polyring('U'+str(self.positives.index(k)+1))

        # creates the new generator
        newgentuple = tuple(-1 if (gen1[i]== -1) else gen2[gen1[i]] for i in range(len(gen1)))
        newgen = Generator(gena.parent,list(newgentuple))

        return (coefficient, newgen)

    # Sage's built-in algebra constructor wants matrices
    def matrix_representation(self, index):
        g = self.gens[index]

        images = {}
        for i in range(len(self.gens)):
            product = self.multiply(self.gens[i], g)
            if product[0] == 0:
                continue
            images[(i, self.gens.index(product[1]))] = product[0]
        return matrix(self.polyring, len(self.gens), len(self.gens), images)
    
    # for easy testing
    def multiply_in_algebra(self, gen1, gen2):
        g1 = self.algebra.gen(self.gens.index(gen1))
        g2 = self.algebra.gen(self.gens.index(gen2))
        return g1*g2
            
