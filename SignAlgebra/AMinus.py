import copy
from SignAlgebra.Z2PolynomialRing import *
from SignAlgebra.AMinusGen import *


class AMinus:
    def __init__(self, sign_sequence):
        # store the overarching sign sequence, the list of positive indices, and the polynomial ring
        self.ss = sign_sequence
        self.positives = list(i for i in range(len(sign_sequence)) if sign_sequence[i] == 1)
        self.polyring = self.initpolyring()

    # construct the polynomial ring based on the number of positives
    def initpolyring(self):
        np = len(self.positives)
        return Z2PolynomialRing(['U%s' % p for p in range(1, np + 1)])

    def multiply(self, elta, eltb):
        # if you call multiply on two generators, multiply
        agen = isinstance(elta, AMinusGen)
        bgen = isinstance(eltb, AMinusGen)
        assert agen == bgen
        if agen and bgen:
            return self.gen_mult(elta, eltb)
        # else extend linearly
        else:
            eltout = {}
            for keya in elta.keys():
                for keyb in eltb.keys():
                    genprod = self.gen_mult(keya, keyb)
                    eltout = self.add(eltout, genprod)
            return AMinus.clean(eltout)

    def diff(self, elt):
        # if a generator, return the differential of that generator
        if isinstance(elt, AMinusGen):
            return AMinus.clean(self.gen_diff(elt))
        # else sum the differentials of all the strand diagrams
        else:
            eltout = {}
            for key in elt.keys():
                diff = self.gen_diff(key)
                diff = AMinus.scale(diff, elt[key])
                eltout = self.add(eltout, diff)
            return AMinus.clean(eltout)

    # scale an arbitrary element by some coefficient
    @staticmethod
    def scale(elt, coeff):
        # assert coeff in self.polyring
        if isinstance(elt, AMinusGen):
            return AMinus.clean({elt: coeff})
        else:
            for key in elt.keys():
                elt[key] = elt[key] * coeff
            return AMinus.clean(elt)

    def add(self, elta, eltb):
        # cases. if a or b is a generator, make it an arbitrary element with coeff 1 then add
        agen = isinstance(elta, AMinusGen)
        bgen = isinstance(eltb, AMinusGen)
        if agen and bgen:
            return self.add({elta: 1}, {eltb: 1})
        elif agen:
            return self.add({elta: 1}, eltb)
        elif bgen:
            return self.add(elta, {eltb: 1})
        # add two elements. create output, add values that are common to both, values just in a, and values just in b
        else:
            eltout = {}
            bkeys = eltb.keys()
            for key in elta.keys():
                if key in bkeys:
                    eltout[key] = elta[key] + eltb[key]
                    bkeys.remove(key)
                else:
                    eltout[key] = elta[key]

            for key in bkeys:
                eltout[key] = eltb[key]

            return AMinus.clean(eltout)

    def gen_mult(self, gena, genb):
        gena = dict(gena.strands)
        genb = dict(genb.strands)
        # check if the ends don't match
        if set(gena.values()) != set(genb.keys()):
            return {}

        # check if black strands double cross
        for keya1 in gena.keys():
            for keya2 in gena.keys():
                if (genb[gena[keya1]] > genb[gena[keya2]] and gena[keya1] < gena[keya2]) \
                        or (genb[gena[keya1]] < genb[gena[keya2]] and gena[keya1] > gena[keya2]):
                    return {}

        # count double crossed orange strands
        coeff = self.polyring.one()
        for key in gena.keys():
            if gena[key] > key:
                checkrange = range(max(i, genb[gena[key]]), gena[key])
            else:
                checkrange = range(gena[key], min(key, genb[gena[key]]))

            for i in checkrange:
                if i not in self.positives:
                    return {}
                else:
                    coeff = coeff * self.polyring['U' + str(self.positives.index(i) + 1)]

        # construct the new generator
        newgen = {}
        for key in gena.keys():
            newgen[key] = genb[gena[key]]

        return {AMinusGen(newgen): coeff}

    def gen_diff(self, gen):
        gen = dict(gen.strands)
        if gen is None:
            return {AMinusGen({}): 0}
        # find all strands that cross, and resolve them
        else:
            eltout = {}
            keys = gen.keys()
            for keya in keys:
                for keyb in keys:
                    if keyb < keya and gen[keyb] > gen[keya]:
                        eltout = self.add(eltout, self.resolve(gen, keya, keyb))
            return eltout

    def resolve(self, gen, i, j):
        keys = gen.keys()
        # check if black strands double cross
        for key in set(keys) & set(range(j, i)):
            if gen[i] < gen[key] < gen[j]:
                return {AMinusGen({}): 0}

        # construct the output generator
        output = copy.deepcopy(gen)
        output[i] = gen[j]
        output[j] = gen[i]
        output = AMinusGen(output)

        # calculate the appropriate coefficient
        checkrange = range(max(gen[i], j), min(i, gen[j]))
        coeff = self.polyring.one()
        for i in checkrange:
            if i not in self.positives:
                return {output: 0}
            else:
                coeff = coeff * self.polyring['U' + str(self.positives.index(i) + 1)]

        return {output: coeff}

    @staticmethod
    def clean(elt):
        if isinstance(elt, AMinusGen):
            return elt
        else:
            for key in list(elt.keys()):
                if elt[key] == 0:
                    del elt[key]
            return elt

    def twoalexander(self, elt):
        if isinstance(elt, AMinusGen):
            return self.twoalex_gen(elt, self.polyring.one())
        else:
            firstkey = elt.keys()[0]
            twoalex = self.twoalex_gen(firstkey, elt[firstkey])
            for key in elt.keys():
                if self.twoalex_gen(key, elt[key]) != twoalex:
                    return None
            return twoalex

    def twoalex_gen(self, gen, coeff):
        gen = dict(gen.strands)
        if gen is None:
            return None
        twoalex = -2 * coeff.degree()
        for key in gen.keys():
            if key < gen[key]:
                checkrange = range(key, gen[key])
            else:
                checkrange = range(gen[key], key)
            for k in checkrange:
                twoalex = twoalex + (-1) * self.ss[k]

        return twoalex

    def maslov(self, elt):
        if isinstance(elt, AMinusGen):
            return self.maslov_gen(elt, self.polyring.one())
        else:
            firstkey = elt.keys()[0]
            maslov = self.maslov_gen(firstkey, elt[firstkey])
            for key in elt.keys():
                if self.maslov_gen(key, elt[key]) != maslov:
                    return None
            return maslov

    def maslov_gen(self, gen, coeff):
        gen = dict(gen.strands)
        if gen is None:
            return None
        maslov = -2 * coeff.degree()
        for key in gen.keys():
            for keyb in gen.keys():
                if keyb < key and gen[keyb] > gen[key]:
                    maslov = maslov + 1
            if key < gen[key]:
                checkrange = range(key, gen[key])
            else:
                checkrange = range(gen[key], key)
            for k in checkrange:
                if k in self.positives:
                    maslov = maslov - 1

        return maslov
