import numpy as np
import density, estimators, rates, helper

def divergences_rate(est_type, ns, ss, d=1, iters=50, fast=True):
    """
    Experiment for plotting rate of convergence of divergence estimator on Renyi and Tsallis divergence.
    Input:
    est_type -- which estimator should we use. "linear" "plugin" and "quadratic" are supported.
    ns -- a list of sample sizes to evaluate on. Try np.logspace(1, 4, 20)
    ss -- a list which smoothness parameters should we evaluate on?
    d -- the dimension
    iters -- how many iterations should we use.
    fast -- should we use fast numeric integration or slow numeric integration

    Returns nothing but writes log files with the results of the experiments to ./data/
    """
    E = None
    if est_type == "plugin":
        E = estimators.PluginEstimator
    elif est_type == "linear":
        E = estimators.LinearEstimator
    elif est_type == "quadratic":
        E = estimators.QuadraticEstimator
    else:
        print "Estimator %s not supported" % (est_type)
        return

    for s in ss:
        print "s = %s" % (str(s))
        if d == 1:
            Dp = density.UniTrigDensity(s, 1)
            Dq = density.UniTrigDensity(s, 1)
        else:
            Dp = density.TrigDensity(s, 1, d)
            Dq = density.TrigDensity(s, 1, d)
        (new_ns, ms, vs) = rates.renyi_rate(E, Dp, Dq, ns, alpha=0.5, iters=iters, fast=fast)
        f = open("./data/%s_renyi_error_d=%d_s=%s.out" % (est_type, d, str(s)), "w")
        f.write("ns " + " ".join([str(n) for n in ns]) + "\n")
        f.write("ms " + " ".join([str(m) for m in ms]) + "\n")
        f.write("vs " + " ".join([str(v) for v in vs]))
        f.close()
        (new_ns, ms, vs) = rates.tsallis_rate(E, Dp, Dq, ns, alpha=0.5, iters=iters, fast=fast)
        f = open("./data/%s_tsallis_error_d=%d_s=%s.out" % (est_type, d, str(s)), "w")
        f.write("ns " + " ".join([str(n) for n in ns]) + "\n")
        f.write("ms " + " ".join([str(m) for m in ms]) + "\n")
        f.write("vs " + " ".join([str(v) for v in vs]))
        f.close()
    return

def l2_rate_experiment(ns, ss, d=1, iters=50, fast=True):
    """
    Experiment for plotting rate of convergence of divergence estimator on L_2^2 divergence.
    Input:
    est_type -- which estimator should we use. "linear" "plugin" and "quadratic" are supported.
    ns -- a list of sample sizes to evaluate on. Try np.logspace(1, 4, 20)
    ss -- a list which smoothness parameters should we evaluate on?
    d -- the dimension
    iters -- how many iterations should we use.
    fast -- should we use fast numeric integration or slow numeric integration

    Returns nothing but writes log files with the results of the experiments to ./data/
    """
    for s in ss:
        print "s = %s" % (str(s))
        if d == 1:
            Dp = density.UniTrigDensity(s, 1)
            Dq = density.UniTrigDensity(s, 1)
        else:
            Dp = density.TrigDensity(s, 1, d)
            Dq = density.TrigDensity(s, 1, d)
        (new_ns, ms, vs) = rates.l2_rate(Dp, Dq, ns, iters=iters, fast=fast)
        f = open("./data/l2_error_d=%d_s=%s.out" % (d, str(s)), "w")
        f.write("ns " + " ".join([str(n) for n in ns]) + "\n")
        f.write("ms " + " ".join([str(m) for m in ms]) + "\n")
        f.write("vs " + " ".join([str(v) for v in vs]))
        f.close()
    return

if __name__=="__main__":
#     ss = np.arange(0.25, 2.1, 0.25)
    ss = [1.0, 2.0]
    ns = np.logspace(1, 4.0, 30)

    divergences_rate("linear", ns, ss, iters=50)
    l2_rate_experiment(ns, ss)
