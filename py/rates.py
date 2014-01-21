import numpy as np
import density, kde
import matplotlib.pyplot as plt

def kde_rate_1D(D, ns, ps, iters=10):
    """
    Compute the kde_rate of convergence on the density D as a function on n.
    The rate will be computed in ell_p^p norm for each p
    """

    ms = [[] for p in ps]
    vs = [[] for p in ps]
    pts = np.matrix(np.arange(0, 1, 0.001)).T

    for n in ns:
        sub_scores = [[] for p in ps]
        for i in range(iters):
            data = D.sample(n)
            K = kde.KDE(data, D.s)
            [sub_scores[i].append(K.kde_error(pts, D, ps[i])) for i in range(len(ps))]
        print sub_scores
        [ms[i].append(np.mean(sub_scores[i])) for i in range(len(ps))]
        [vs[i].append(np.std(sub_scores[i])) for i in range(len(ps))]
    
    return (ns, ms, vs)
