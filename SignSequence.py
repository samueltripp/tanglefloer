class SignSequence:
    def __init__(self, sequence):
        for i in sequence:
            assert (i == 1) or (i == -1)
        self.sequence = sequence
        self.degree = len(sequence)

    def degree(self):
        return self.degree

    def sequence(self):
        return self.sequence
