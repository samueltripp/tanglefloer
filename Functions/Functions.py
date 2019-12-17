import itertools


# helper functions

# returns a list of all injections, modeled as dicts
def injections(source, target):
    if len(source) == 0:
        return [{}]
    if len(target) == 0:
        return []
    output = []
    for index, image in enumerate(target):
        output.extend([copy_and_add(injection, source[0], image) for injection
                       in injections(source[1:], target[:index] + target[index + 1:])])
    return output


# returns a list of all partial bijections, modeled as dicts
def partial_bijections(source, target):
    output = []
    for sublist in sublists(source):
        output.extend(injections(sublist, target))
    return output


def invert_injection(inj):
    return {t: s for s, t in inj.items()}


# returns a list of all sublists of the given list
def sublists(a_list):
    output = []
    for i in range(len(a_list) + 1):
        output.extend(itertools.combinations(a_list, i))
    return output


# creates a copy of the given dictionary d with (k,v) added
def copy_and_add(d, k, v):
    new_dict = dict(d)
    new_dict[k] = v
    return new_dict


# creates a copy of the given dictionary d with (-,v) removed
def copy_and_remove(d, value_to_remove):
    return {k: v for k, v in d.items() if v != value_to_remove}


# creates a point with flipped y-orientation
def invert(x, y):
    return x, -y


# an alternate __str__() for dictionaries that sorts the keys first
def dict_to_sorted_string(d) -> str:
    if len(d) == 0:
        return '{}'
    out = '{'
    for k in sorted(d.keys()):
        out += str(k) + ': ' + str(d[k]) + ', '
    return out[:-2] + '}'
