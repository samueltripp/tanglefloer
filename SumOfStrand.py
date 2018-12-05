class SumOfStrand:
    def __init__(self,dict):
        self.dict = dict

    def dict(self):
        return dict

    def add(self,sd,coeff):
        dict = self.dict
        if sd in dict:
            dict[sd] = dict[sd]+coeff
        else :
            dict[sd] = coeff

    def pop(self,key):
        self.dict.pop(key,None)
