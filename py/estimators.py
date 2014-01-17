import numpy as np
import kde, density
import itertools

class PluginEstimator(object):
    """
    Plugin estimator for functionals of the form \int
    p^\alpha(x)q^\beta(x). We construct kernel density estimators for
    p,q$ and then perform numerically integrate to evalute \int \phat^\alpha(x) \qhat^\beta(x)
    """
    def __init__(self, pdata, qdata, alpha, beta, s):
        self.Kp = kde.KDE(pdata, s)
        self.Kq = kde.KDE(qdata, s)
        self.alpha = alpha
        self.beta = beta
        self.s = s
        
    def eval(self, n):
        """
        Evaluate the estimator by performing numeric integration on the n^d grid of points.
        """
        nr = range(n)
        vals = []
        for i in itertools.combinations_with_replacement(nr, self.Kp.d):
            pt = np.matrix(i)/float(n)
            vals = self.Kp.eval(pt)[0]**self.alpha * self.Kq.eval(pt)[0]**self.beta
        return np.mean(vals)

class LinearEstimator(object):
    pass

class QuadraticEstimator(object):
    pass

class Truth(object):
    def __init__(self, p, q, alpha, beta):
        self.p = p
        self.q = q
        self.alpha = alpha
        self.beta = beta

    def eval(self, n):
        nr = range(n)
        vals = []
        for i in itertools.combinations_with_replacement(nr, self.p.d):
            pt = np.matrix(i)/float(n)
            vals = self.p.eval(pt)[0]**self.alpha * self.q.eval(pt)[0]**self.beta
        return np.mean(vals)
        


if __name__=='__main__':
    print "Generating two uni-trig-densities"
    Dp = density.UniTrigDensity(2, 2)
    Dq = density.UniTrigDensity(2, 2)    
    alpha = 1
    beta = 1

    T = Truth(Dp, Dq, alpha, beta)
    print "Numeric integration of true functional"
    tr = T.eval(100)
    print "T(p,q) approx %0.2f" % (tr)

    scores = []
    for n in range(100, 1000, 100):
        print "n = %d" % (n)
        sub_scores = []
        for i in range(10):
            pdata = Dp.sample(n)
            qdata = Dq.sample(n)            
            PL = PluginEstimator(pdata, qdata, alpha, beta, 2)
            sc = PL.eval(100)
            print "Score: %0.2f" % (sc)
            sub_scores.append(sc)
        scores.append(sub_scores)
    
    print scores
