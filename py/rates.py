import numpy as np
from helper import *
import density, kde, estimators
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


def estimator_rate(Est, Dp, Dq, ns, alpha=1, beta=1, iters=10):
    """
    Compute the rate of convergenc of the estimator Est on the densities Dp, Dq
    The rate will be computed in \ell_1 norm.
    """
    ms = []
    vs = []
    T = estimators.Truth(Dp, Dq, alpha, beta)
    truth = T.eval()

    for n in ns:
        sub_scores = []
        for i in range(iters):
            pdata = Dp.sample(n)
            qdata = Dq.sample(n)
            E = Est(pdata, qdata, alpha, beta, Dp.s)
            val = E.eval()
            sub_scores.append(np.abs(val-truth))

        ms.append(np.mean(sub_scores))
        vs.append(np.std(sub_scores))
        print "n = %d, av_er = %0.2f, std_er = %0.4f" % (n, np.mean(sub_scores), np.std(sub_scores))
    return (ns, ms, vs)


# To rescale a rate of n^{-\gamma}, plot ns versus 
# [np.exp(np.log(ms[i]) + gamma*np.log(ns[i])) for i in range(len(ns))]
# And you should see the constant \gamma come out.

# To empirically verify the parameter \gamma, plot ns versus
# [- np.log(ms[i])/np.log(ns[i]) for i in range(len(ns))]
# And you should see the line approach \gamma. 
