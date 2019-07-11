# this only exists because I want to use strand diagrams as keys in algebra elements. strand diagrams are dictionaries, however, and you can't use dictionaries as keys in dictionaries. hence, wrap them in a class. 

class AMinusGen:
    def __init__(self,dict):
        self.dict = dict

    def __eq__(self,other):
        return self.dict == other.dict
    
    def __hash__(self):
        return hash(str(self))
