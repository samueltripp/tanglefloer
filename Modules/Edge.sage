# represents the action of some delta_1 on some pair of generators
class Edge:
    # source, target - elements of gens
    # a_coefficient - element of A
    # m_cofficient - element of k
    # b_coefficient - element of B
    # 
    # condition: delta_1(source) = a_coefficient (X) (m_coefficient * target) (X) b_coefficient
    def __init__(self, source, target, a_coefficient, b_coefficient, m_coefficient):
        self.source = source
        self.target = target
        self.a_coefficient = a_coefficient
        self.b_coefficient = b_coefficient
        self.m_coefficient = m_coefficient
