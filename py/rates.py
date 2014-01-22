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

    for n in ns:
        sub_scores = [[] for p in ps]
        for i in range(iters):
            data = D.sample(n)
            K = kde.KDE(data, D.s)
            [sub_scores[i].append(K.kde_error(D, ps[i])) for i in range(len(ps))]
        print "n = %d " % n + " ".join([str(ps[i]) + " = %0.2f" % np.mean(sub_scores[i]) for i in range(len(ps))])
        [ms[i].append(np.mean(sub_scores[i])) for i in range(len(ps))]
        [vs[i].append(np.std(sub_scores[i])) for i in range(len(ps))]
    
    return (ns, ms, vs)


def estimator_rate(Est, Dp, Dq, ns, alpha=1, beta=1, iters=10, fast=True):
    """
    Compute the rate of convergenc of the estimator Est on the densities Dp, Dq
    The rate will be computed in \ell_1 norm.
    """
    ms = []
    vs = []
    T = estimators.Truth(Dp, Dq, alpha, beta)
    truth = T.eval(fast=fast)

    for n in ns:
        sub_scores = []
        for i in range(iters):
            pdata = Dp.sample(n)
            qdata = Dq.sample(n)
            E = Est(pdata, qdata, alpha, beta, Dp.s)
            val = E.eval(fast=fast)
            sub_scores.append(np.abs(val-truth))

        ms.append(np.mean(sub_scores))
        vs.append(np.std(sub_scores))
        print "n = %d, av_er = %0.2f, std_er = %0.4f" % (n, np.mean(sub_scores), np.std(sub_scores))
    return (ns, ms, vs)


def plot_estimator_rate(ns, ms, vs, gamma):
    """
    Generate three plots. 
    One is an error bar plot of the error of the estimator.
    The next is n vs \exp\{\log(err(n)) + \gamma \log n\} which should pull out the constant in the rate.
    The third is n vs -\log(err(n)/\log(n) which should approach \gamma
    """
    fig = plt.figure(figsize=(15, 5))
    ax1 = fig.add_subplot(131)
    ax1.errorbar(ns, ms, vs)
    ax1.set_xlabel("Number of samples (n)")
    ax1.set_ylabel("Error |That - T|")
    ax2 = fig.add_subplot(132)
    ax2.plot(ns, [ms[i]*ns[i]**gamma for i in range(len(ns))])
    ax2.set_xlabel("Number of samples (n)")
    ax2.set_ylabel("Error*n^{\gamma}")
    ax3 = fig.add_subplot(133)
    ax3.plot(ns, [-np.log(ms[i])/np.log(ns[i]) for i in range(len(ns))])
    ax3.set_xlabel("Number of samples (n)")
    ax3.set_ylabel("-log(error)/log(n)")
    plt.show()

# To rescale a rate of n^{-\gamma}, plot ns versus 
# [np.exp(np.log(ms[i]) + gamma*np.log(ns[i])) for i in range(len(ns))]
# And you should see the constant \gamma come out.

# To empirically verify the parameter \gamma, plot ns versus
# [- np.log(ms[i])/np.log(ns[i]) for i in range(len(ns))]
# And you should see the line approach \gamma. 
