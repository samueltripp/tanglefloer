

# this only exists because I want to use strand diagrams as keys in algebra elements.
# strand diagrams are dictionaries, however, and you can't use dictionaries as keys in dictionaries.
# hence, wrap them in a class.
class AMinusGen:
    def __init__(self, strands):
        self.strands = frozenset(strands.items())

    def __eq__(self, other):
        return self.strands == other.strands
    
    def __hash__(self):
        return hash(str(self.strands))
