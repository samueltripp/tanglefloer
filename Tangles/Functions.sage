import itertools


# helper functions


# returns a list of all injections, modeled as dicts
def injections(source, target):
	if len(source) == 0:
		return [{}]
	elif len(target) < len(source):
		return []
	output = []
	for index,image in enumerate(target):
		for inj in injections(source[1:], target[:index] + target[index + 1:]):
			inj[source[0]] = image
			output.append(inj)
	return output


# returns a list of all partial bijections, modeled as dicts
def partial_bijections(source, target):
	output = []
	for sublist in sublists(source, min(len(source), len(target))):
		output.extend(injections(sublist, target))
	return output


# returns a list of all sublists of the given list of size at most k
def sublists(a_list, k):
	output = []
	for i in range(k + 1):
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


# # creates a point with flipped y-orientation
# def c(x,y):
#     return (x, -y)
