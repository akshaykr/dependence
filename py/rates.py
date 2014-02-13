import density
import estimators
import kde
import numpy as np

from helper import numeric_integration
from plotting import (
    plot_estimator_rate,
    plot_log_log,
)
from tests import test_linear_functional

def kde_rate(D, ns, ps, iters=10):
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
    truth = T.eval(fast=True)
#     truth = 0.25 * (2 + 5.0* np.arctan(1.0/2))

    for n in ns:
        sub_scores = []
        for i in range(iters):
            pdata = Dp.sample(n)
            qdata = Dq.sample(n)
            E = Est(pdata, qdata, alpha, beta, Dp.s)
            val = E.eval(fast=fast)
            sub_scores.append(np.abs(val-truth))

        sub_scores.sort();
        sub_scores = sub_scores[int(0.2*iters): int(0.8*iters)];
        ms.append(np.mean(sub_scores))
        vs.append(np.std(sub_scores))
        print "n = %d, av_er = %0.2f, std_er = %0.4f" % (n, np.mean(sub_scores), np.std(sub_scores))
    return (ns, ms, vs)

def renyi_rate(Est, Dp, Dq, ns, alpha=0.5, iters=50, fast=True):
    """
    Rate of convergence experiment for Renyi-alpha divergence.
    """
    ms = []
    vs = []
    T = estimators.Truth(Dp, Dq, alpha, 1-alpha)
    truth = 1.0/(alpha-1) * np.log(T.eval(fast=False))

    for n in ns:
        sub_scores = []
        for i in range(iters):
            pdata = Dp.sample(n)
            qdata = Dq.sample(n)
            E = Est(pdata, qdata, alpha, 1-alpha, Dp.s)
            val = E.eval(fast=fast)
            val = 1.0/(alpha-1) * np.log(val)
            sub_scores.append(np.abs(val-truth))

        sub_scores.sort();
        sub_scores = sub_scores[int(0.2*iters): int(0.8*iters)];
        ms.append(np.mean(sub_scores))
        vs.append(np.std(sub_scores))
        print "n = %d, av_er = %0.2f, std_er = %0.4f" % (n, np.mean(sub_scores), np.std(sub_scores))
    return (ns, ms, vs)

def tsallis_rate(Est, Dp, Dq, ns, alpha=0.5, iters=50, fast=True):
    """
    Rate of convergence experiment for tsallis-alpha divergence.
    """
    ms = []
    vs = []
    T = estimators.Truth(Dp, Dq, alpha, 1-alpha)
    truth = 1.0/(alpha-1) * (T.eval(fast=False) - 1)

    for n in ns:
        sub_scores = []
        for i in range(iters):
            pdata = Dp.sample(n)
            qdata = Dq.sample(n)
            E = Est(pdata, qdata, alpha, 1-alpha, Dp.s)
            val = E.eval(fast=fast)
            val = 1.0/(alpha-1) * (val - 1)
            sub_scores.append(np.abs(val-truth))

        sub_scores.sort();
        sub_scores = sub_scores[int(0.2*iters): int(0.8*iters)];
        ms.append(np.mean(sub_scores))
        vs.append(np.std(sub_scores))
        print "n = %d, av_er = %0.2f, std_er = %0.4f" % (n, np.mean(sub_scores), np.std(sub_scores))
    return (ns, ms, vs)

def l2_rate(Dp, Dq, ns, iters=50, fast=True):
    """
    Collect data for L_2^2 divergence estimation of the distributions Dp, Dq over the set of n-values.
    Repeat trials iters times.
    """
    ms = []
    vs = []
    x = np.matrix(np.arange(0, 1, 0.001)).T
    truth = np.mean(np.array((Dp.eval(x)- Dq.eval(x)))**2)

    truth = numeric_integration(lambda x:
                                (Dp.eval(np.matrix(x))- Dq.eval(np.matrix(x)))**2,
                                [0], [1])
    for n in ns:
        sub_scores = []
        for i in range(iters):
            pdata = Dp.sample(n)
            qdata = Dq.sample(n)
            E = estimators.L2Estimator(pdata, qdata, Dp.s)
            val = E.eval(fast=fast)
            sub_scores.append(np.abs(val-truth))

        sub_scores.sort();
        sub_scores = sub_scores[int(0.2*iters): int(0.8*iters)];
        ms.append(np.mean(sub_scores))
        vs.append(np.std(sub_scores))
        print "n = %d, av_er = %0.2f, std_er = %0.4f" % (n, np.mean(sub_scores), np.std(sub_scores))
    return (ns, ms, vs)

if __name__ == '__main__' :
  D = density.UniTrigDensity(2,2);
  (ns, ms, vs) = test_linear_functional(D, np.logspace(1, 4, 20), iters = 50);
  plot_estimator_rate(ns, ms, vs, 0.5)
  (m, b) = plot_log_log(ns, ms)
  print m, b

# To rescale a rate of n^{-\gamma}, plot ns versus
# [np.exp(np.log(ms[i]) + gamma*np.log(ns[i])) for i in range(len(ns))]
# And you should see the constant \gamma come out.

# To empirically verify the parameter \gamma, plot ns versus
# [- np.log(ms[i])/np.log(ns[i]) for i in range(len(ns))]
# And you should see the line approach \gamma.
