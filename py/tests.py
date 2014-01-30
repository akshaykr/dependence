import numpy as np
import estimators
from helper import *

def test_integration(Dp, Dq, ns, alpha=0.5):
    """
    Test the accuracy of the integrator for Dp, Dq.
    Use n rectangles for each value of ns
    """
    T = estimators.Truth(Dp, Dq, alpha, 1-alpha)
    truth = T.eval(fast=False)

    ms = []
    for n in ns:
        val = T.eval(fast=True, pts=n)
        ms.append(np.abs(truth-val))
    return (ns, ms)

def test_quadratic_term_estimator(Dp, ns, iters=10, fast=True):
    """
    Test the quadratic term estimator.
    """
    ms = []
    vs = []
    T = fast_integration(lambda x: np.array(Dp.eval(x))**2, [0], [1])
    print "Truth: %f" % (T)
    for n in ns:
        sub_scores = []
        for i in range(iters):
            pdata = Dp.sample(n)
            Q = estimators.QuadraticEstimator(pdata, pdata, 0.5, 0.5, Dp.s)
            val = Q.quad_term_slow(lambda x: 1, pdata)
            val2 = Q.quad_term_fast(lambda x: 1, pdata)
            print "truth = %0.3f, fast = %0.3f, slow = %0.3f" % (T, val2, val)
            sub_scores.append(np.abs(val-T))
        sub_scores.sort();   
        sub_scores = sub_scores[int(0.2*iters): int(0.8*iters)];
        ms.append(np.mean(sub_scores))
        vs.append(np.std(sub_scores))
        print "n = %d, av_er = %0.2f, std_er = %0.4f" % (n, np.mean(sub_scores), np.std(sub_scores))
    return (ns, ms, vs)

def test_linear_functional(Dp, ns, iters=10, fast=True):
    """
    Test the linear functional estimator
    """
    ms = [];
    vs = [];
    T = numeric_integration(lambda x: np.array(Dp.eval(x)) * np.array(x), [0], [1]); 
    print "True Expectation: %f\n" % (T)
    for n in ns:
        sub_scores = [];
        for i in range(iters):
            pdata = Dp.sample(n);
            val = np.mean(pdata, axis=0)

            sub_scores.append(np.abs(val-T))
        sub_scores.sort();
        sub_scores = sub_scores[int(0.2*iters): int(0.8*iters)];
        ms.append(np.mean(sub_scores))
        vs.append(np.std(sub_scores))
        print "n = %d, av_er = %0.2f, std_er = %0.4f" % (n, np.mean(sub_scores), np.std(sub_scores))
    return (ns, ms, vs)

def test_linear_functional2(Dp, Dq, ns, alpha=0.5, beta=0.5, iters=10, fast=True):
    """
    Another test case for the linear functional
    """
    ms = []
    vs = []
    T = numeric_integration(lambda x: np.multiply(np.power(Dp.eval(x), alpha), np.power(Dq.eval(x), beta)), [0], [1])
    print "True Expectation: %f\n" % (T)
    for n in ns:
        sub_scores = [];
        for i in range(iters):
            pdata = Dp.sample(n);
            qdata = Dq.sample(n)
            fn1 = lambda x: np.multiply(np.power(Dp.eval(x), alpha-1), np.power(Dq.eval(x), beta))
            val1 = alpha* np.mean(fn1(pdata))
            fn2 = lambda x: np.multiply(np.power(Dp.eval(x), alpha), np.power(Dq.eval(x), beta-1))
            val2 = beta* np.mean(fn1(qdata))
            sub_scores.append(np.abs(val1+val2-T))
        sub_scores.sort();
        sub_scores = sub_scores[int(0.2*iters): int(0.8*iters)];
        ms.append(np.mean(sub_scores))
        vs.append(np.std(sub_scores))
        print "n = %d, av_er = %0.2f, std_er = %0.4f" % (n, np.mean(sub_scores), np.std(sub_scores))
    return (ns, ms, vs)

def test_bilinear_term_estimator(Dp, Dq, ns, iters=10, fast=True):
    """
    Test the bilinear functional estimator
    """
    ms = []
    vs = []
    T = fast_integration(lambda x: np.array(Dp.eval(x))*np.array(Dq.eval(x)), [0], [1])
    print "Truth: %f" % (T)
    for n in ns:
        sub_scores = []
        for i in range(iters):
            pdata = Dp.sample(n)
            qdata = Dq.sample(n)
            Q = estimators.QuadraticEstimator(pdata, pdata, 0.5, 0.5, Dp.s)
            val = Q.bilinear_term_slow(lambda x: 1, pdata, qdata)
            val2 = Q.bilinear_term_fast(lambda x: np.matrix(np.ones((x.shape[0], 1))), pdata, qdata)
            print "truth = %0.3f, slow = %0.3f, fast = %0.3f" % (T, val2, val)
            sub_scores.append(np.abs(val-T))
        sub_scores.sort();   
        sub_scores = sub_scores[int(0.2*iters): int(0.8*iters)];
        ms.append(np.mean(sub_scores))
        vs.append(np.std(sub_scores))
        print "n = %d, av_er = %0.2f, std_er = %0.4f" % (n, np.mean(sub_scores), np.std(sub_scores))
    return (ns, ms, vs)

