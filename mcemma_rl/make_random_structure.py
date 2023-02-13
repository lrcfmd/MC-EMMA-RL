from .all import *
import sys
from random import choice
from random import shuffle
import numpy
import time


def get_composition(AMods, BMods, struct):
    l = len(struct)

    unique_elements = []
    for i in range(len(AMods)):
        atoms = AMods[i]
        for j in range(len(atoms)):
            if atoms[j].symbol not in unique_elements:
                unique_elements.append(atoms[j].symbol)

    for i in range(len(BMods)):
        atoms = BMods[i]
        for j in range(len(atoms)):
            if atoms[j].symbol not in unique_elements:
                unique_elements.append(atoms[j].symbol)

    # test to see if the module set has the target composition
    symbols = {}
    for i in unique_elements:
        symbols[i] = 0.

    for i in range(l):  # go through and add the composition into symbols
        atoms = AMods[struct[i][0]].copy()
        syms = atoms.get_chemical_symbols()
        for j in range(len(syms)):
            symbols[syms[j]] += 1

        atoms = BMods[struct[i][1]].copy()
        syms = atoms.get_chemical_symbols()
        for j in range(len(syms)):
            symbols[syms[j]] += 1

    # print(symbols)
    # check to see if any elements have a value of zero
    keys = list(symbols.keys())

    # find element with smallest value
    nums = []
    for i in range(len(keys)):
        # if symbols[keys[i]] == 0:
        #     struct = False
        #     return struct  # attempt failed, return to main code
        nums.append(symbols[keys[i]])

    smallest = keys[nums.index(min([n for n in nums if n > 0]))]
    smallest = symbols[smallest]
    # print(smallest)
    new_symbols = {}

    for i in range(len(keys)):
        new_symbols[keys[i]] = float(symbols[keys[i]] / smallest)

    return new_symbols


# replaces layers_n layers in row by random substructure with the same composition
def replace_modules(AMods, BMods, struct, layers_n=1, tmax=100000):
    sl = len(struct)

    # choose the start index for the substructure
    start_index = choice(range(sl))

    attempts = 0
    while attempts < sl:

        # choose the list of layes_n indices in row
        # (it can be 4 0 1 for stack length 5 and layers_n=3 and start_index=4)
        indices = list(range(start_index, min(start_index + layers_n, sl)))
        if len(indices) < layers_n:
            indices = indices + list(range(0, layers_n - len(indices)))

        # print("Start index is %d, layers_n is %d, indices are %s" % (start_index, layers_n, str(indices)))

        substruct = struct[indices]
        # print('Struct: ' + str(struct))
        # print('Substruct: ' + str(substruct))

        composition = get_composition(AMods, BMods, substruct)

        # make tmax attempts to create a randum structure for the composition corresponding to the substructure
        z = 0
        while z < tmax/sl:
            new_substruct = make_random_structure(AMods, BMods, composition, [len(indices)])
            if type(new_substruct) != numpy.ndarray:
                z += 1
            elif np.array_equal(substruct, new_substruct):
                new_substruct = False
                z += 1
            else:
                z = tmax/sl

        # if we failed to create a random structure try change start index to try the next possible substructure
        if type(new_substruct) != numpy.ndarray:

            attempts += 1
            print('Attempt %d failed' % attempts)
            # shift start index by one
            if start_index == sl - 1:
                start_index = 0
            else:
                start_index += 1

            continue

        # print('New substruct: ' + str(new_substruct))

        # fill new struct as a hybrid of the old struct and new substruct
        new_struct = struct.copy()
        new_struct[indices] = new_substruct

        # print('New struct: ' + str(new_struct))
        # print('Composition for new struct: ' + str(get_composition(AMods, BMods, new_struct)))

        return new_struct

    print('Replace modules failed with structure ' + str(struct))
    return False


def make_random_structure(AMods, BMods, composition, sl, testing=''):
    ## bin modules by their composition
    if testing == True:
        pass
        # print("\nAMods\n\n",AMods,"\nBMods\n\n",BMods,"\ncomposition\n\n",composition,"\nsl\n\n",sl)

    acomp = {}
    bcomp = {}
    for i in range(len(BMods)):
        try:
            bcomp[BMods[i].get_chemical_formula()].append(i)
        except KeyError:
            bcomp[BMods[i].get_chemical_formula()] = [i]

    for i in range(len(AMods)):
        try:
            acomp[AMods[i].get_chemical_formula()].append(i)
        except KeyError:
            acomp[AMods[i].get_chemical_formula()] = [i]

    param1 = list(acomp.keys())  ### A layer compositions
    param2 = list(bcomp.keys())  ### B layer compositions

    # print("\n\n",acomp)
    # print("\n\n",bcomp)

    if testing == True:
        pass

    # make sure we've got a list of the unique elements in the calculation

    unique_elements = []
    for i in range(len(AMods)):
        atoms = AMods[i]
        for j in range(len(atoms)):
            if atoms[j].symbol not in unique_elements:
                unique_elements.append(atoms[j].symbol)

    for i in range(len(BMods)):
        atoms = BMods[i]
        for j in range(len(atoms)):
            if atoms[j].symbol not in unique_elements:
                unique_elements.append(atoms[j].symbol)

    # collect a module set
    module_set = {'A': [], 'B': []}  # the two module sequences we'll form into a structure
    # print(composition)
    make = False
    l = choice(sl)
    # l=10
    # print(l)
    for i in range(l):
        mod = choice(param1)  # choose A mod to use
        # see how many equivalent versions of it there are
        if len(acomp[mod]) == 1:
            module_set['A'].append(acomp[mod][0])

        if len(acomp[mod]) > 1:
            module_set['A'].append(choice(acomp[mod]))

        mod = choice(param2)  # choose B mod to use
        if len(bcomp[mod]) == 1:
            module_set['B'].append(bcomp[mod][0])

        if len(bcomp[mod]) > 1:
            module_set['B'].append(choice(bcomp[mod]))

    # test to see if the module set has the target composition
    symbols = {}
    for i in unique_elements:
        symbols[i] = 0.

    for i in range(l):  # go through and add the composition into symbols
        atoms = AMods[module_set['A'][i]].copy()
        syms = atoms.get_chemical_symbols()
        for j in range(len(syms)):
            symbols[syms[j]] += 1

        atoms = BMods[module_set['B'][i]].copy()
        syms = atoms.get_chemical_symbols()
        for j in range(len(syms)):
            symbols[syms[j]] += 1

    # print(symbols)
    # check to see if any elements have a value of zero
    keys = list(symbols.keys())

    # find element with smallest value
    nums = []
    for i in range(len(keys)):
        # if symbols[keys[i]] == 0:
        #     struct = False
        #     print(symbols)
        #     return struct  # attempt failed, return to main code
        nums.append(symbols[keys[i]])

    # smallest not null value
    smallest = keys[nums.index(min([n for n in nums if n > 0]))]
    smallest = symbols[smallest]
    # print(smallest)
    new_symbols = {}

    for i in range(len(keys)):
        new_symbols[keys[i]] = float(symbols[keys[i]] / smallest)

    if testing == True:
        print(composition)
        print(new_symbols)
        print(new_symbols == composition)
        sys.exit()
        # print(module_set)

    # print(new_symbols)
    # print(composition)
    # time.sleep(2)
    if new_symbols == composition:
        struct = numpy.zeros([l, 2], dtype=int)
        shuffle(module_set['A'])
        shuffle(module_set['B'])
        for i in range(len(struct)):
            struct[i][0] = module_set['A'][i]
            struct[i][1] = module_set['B'][i]
        return struct  # if we're going to use this method, convert module_set into a struct object then return
    if new_symbols != composition:
        struct = False
        return struct


