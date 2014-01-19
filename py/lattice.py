import numpy as np
import itertools

def lattice(d, limit=10):
    """
    Generator for the integer lattice in d dimensions.
    Not the most efficient thing.
    """
    yield [0 for i in range(d)]
    s = 1
    while True:
        if limit != None and s > limit:
            return
        s += 1
        for i in range(s**d):
            tuple = [i/s**j % s for j in range(0,d)]
            if np.max(tuple) == s-1:
                for g in all_sign_patterns(tuple):
                    yield g
    

def all_sign_patterns(v):
    """
    For a vector v, flip all of the signs of the vector and yield them.
    This is a generator.
    """
    coords = [i for i in range(len(v)) if v[i] != 0]
    for i in range(2**len(coords)):
        s = np.binary_repr(i) + reduce(lambda a,b: a+b, ["0" for k in range(len(coords) - len(np.binary_repr(i)))], "")
        tuple = [int(s[j])*2-1 for j in range(len(s))]
        new_tuple = [0 for i in range(len(v))]
        for i in range(len(coords)):
            new_tuple[coords[i]] = tuple[i]
        yield [v[i]*new_tuple[i] for i in range(len(v))]
