# represents a pair of partial bijections overlaid on an elementary tangle                                                                                
class StrandDiagram:
    def __init__(self, etangle: ETangle, left_strands: Dict, right_strands: Dict):
        self.etangle = etangle
        self.left_strands = left_strands
        self.right_strands = right_strands

    # the idempotent e^D_L                                                                                                                                
    def left_idempotent(self):
        occupied = self.left_strands.keys()
        total = set(range(len(self.etangle.left_algebra.ss)))
        return self.etangle.left_algebra.idempotent(list(total - occupied))

    # the idempotent e^A_R                                                                                                                                
    def right_idempotent(self) -> AMinusElement:
        return self.etangle.right_algebra.idempotent(list(self.right_strands.values()))

    def __repr__(self):
        return str((self.etangle, self.left_strands, self.right_strands))
