class SignSequence:
    def __init__(self, sequence):
        for i in sequence:
            assert (i == 1) or (i == -1)
        self.sequence = sequence
        self.degree = len(sequence)

        # later we need to know, given a positively oriented strand, what count positive strand this is
        temp = 1
        npos = 0
        poscount = [-1]*len(sequence)
        for i in range(0,len(sequence)):
            if sequence[i] == 1:
                poscount[i] = temp
                temp = temp + 1
                npos = npos + 1
        self.poscount = poscount
        self.npos = npos


    def degree(self):
        return self.degree

    def sequence(self):
        return self.sequence

    def poscount(self):
        return self.poscount

    def npos(self):
        return self.npos

    def __str__(self):
        return str(self.sequence)

    def __repr__(self):
        return str(self.sequence)

    def __eq__(self,other):
        if (self.sequence == other.sequence): return True
        else: return False



