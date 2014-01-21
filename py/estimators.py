import numpy as np
import kde, density
import itertools
from helper import *

class PluginEstimator(object):
    """
    Plugin estimator for functionals of the form \int
    p^\alpha(x)q^\beta(x). We construct kernel density estimators for
    p,q$ and then perform numerically integrate to evalute \int \phat^\alpha(x) \qhat^\beta(x)
    """
    def __init__(self, pdata, qdata, alpha, beta, s):
        self.Kp = kde.KDE(pdata, s)
        self.Kq = kde.KDE(qdata, s)
        self.dim = pdata.shape[1]
        self.alpha = alpha
        self.beta = beta
        self.s = s
        
    def eval(self, fast=False):
        """
        Evaluate the estimator by performing numeric integration on the n^d grid of points.
        """
        if self.dim == 1:
            if fast:
                val = fast_integration(lambda x: np.multiply(
                    np.power(self.Kp.eval(np.matrix(x)), self.alpha), 
                    np.power(self.Kq.eval(np.matrix(x)), self.beta)),
                                       [0], [1])
            else:
                val = numeric_integration(lambda x: np.multiply(
                        np.power(self.Kp.eval(np.matrix(x)), self.alpha), 
                        np.power(self.Kq.eval(np.matrix(x)), self.beta)),
                                          [0], [1])
        if self.dim == 2:
            if fast:
                val = fast_integration(lambda x,y: np.multiply(
                        np.power(self.Kp.eval(np.concatenate((np.matrix(x), np.matrix(y)),
                                                             1)), self.alpha),
                        np.power(self.Kq.eval(np.concatenate((np.matrix(x), np.matrix(y)),
                                                             1)), self.beta)),
                                          [0,0], [1,1])
            else:
                val = numeric_integration(lambda x,y: np.multiply(
                        np.power(self.Kp.eval(np.concatenate((np.matrix(x), np.matrix(y)),
                                                             1)), self.alpha),
                        np.power(self.Kq.eval(np.concatenate((np.matrix(x), np.matrix(y)),
                                                             1)), self.beta)),
                                          [0,0], [1,1])
        return val


class LinearEstimator(PluginEstimator):

    def __init__(self, pdata, qdata, alpha, beta, s):
        # Shuffle the data and split it into 2 for density estimation and
        # estimating the other parts.
        np.random.shuffle(pdata);
        np.random.shuffle(qdata);
        num_pdata = pdata.shape[0];
        num_qdata = qdata.shape[0];
        # Call the previous constructor
        PluginEstimator.__init__(self, pdata[0:num_pdata/2, :],
          qdata[0:num_qdata/2, :], alpha, beta, s);
        # Store the 2 splits
        self.p_den_data = pdata[0: num_pdata/2, :];
        self.q_den_data = qdata[0: num_qdata/2, :];
        self.p_est_data = pdata[num_pdata/2 + 1: num_pdata, :];
        self.q_est_data = qdata[num_qdata/2 + 1: num_qdata, :];

    def eval(self):
        # Plugin estimator
        plugin_est = super(LinearEstimator, self).eval();

        # C1 = - (alpha + beta) * int {p0(x)^alpha q0(x)^beta} dx
        if self.dim == 1:
          C1_fn_handle = lambda x: np.multiply( 
            np.power(self.Kp.eval(np.matrix(x)), self.alpha),
            np.power(self.Kq.eval(np.matrix(x)), self.beta) );
          l_limit = [0];
          u_limit = [1];

        elif self.dim == 2:
          C1_fn_handle = lambda x,y : np.multiply( 
            np.power(self.Kp.eval( np.concatenate((np.matrix(x), np.matrix(y)),
                                   1) ), self.alpha),
            np.power(self.Kq.eval( np.concatenate((np.matrix(x), np.matrix(y)),
                                   1) ), self.beta) );
          l_limit = [0, 0];
          u_limit = [1, 1];
        C1 = -(self.alpha + self.beta) * numeric_integration(
                                           C1_fn_handle, l_limit, u_limit);

        # theta^p_{1,1} = alpha * E_{x~p}[p0(x)^(alpha-1) * q0(x)^beta]
        theta_p_11 = self.alpha * np.mean(
            np.multiply(np.power(self.Kp.eval(self.p_est_data), self.alpha - 1),
                        np.power(self.Kq.eval(self.p_est_data), self.beta)
                       ) );

        # theta^q_{1,1} = beta * E_{x~q}[p0(x)^alpha * q0(x)^{beta-1}]
        theta_q_11 = self.beta * np.mean(
            np.multiply(np.power(self.Kp.eval(self.q_est_data), self.alpha),
                        np.power(self.Kq.eval(self.q_est_data), self.beta - 1)
                       ) );
        # Return plugin + C1 + theta_p_11 + theta_q_11
        return plugin_est + theta_p_11 + theta_q_11 + C1;


class QuadraticEstimator(object):
    pass

class Truth(object):
    def __init__(self, p, q, alpha, beta):
        self.p = p
        self.q = q
        self.alpha = alpha
        self.beta = beta

    def eval(self):
        """
        Evaluate the true divergence by performing numeric integration on the n^d grid of points.
        """
        if self.dim == 1:
            val = numeric_integration(lambda x: np.multiply(
                    np.power(self.p.eval(np.matrix(x)), self.alpha), 
                    np.power(self.q.eval(np.matrix(x)), self.beta)),
                    [0], [1])
        if self.dim == 2:
            val = numeric_integration(lambda x,y: np.multiply(
                    np.power(self.p.eval(np.concatenate((np.matrix(x), np.matrix(y)),
                                                          1)), self.alpha),
                    np.power(self.q.eval(np.concatenate((np.matrix(x), np.matrix(y)),
                                                          1)), self.beta)),
                                      [0,0], [1,1])
        return val
        


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

    pl_scores = []
    lin_scores = [];
    for n in range(100, 1000, 100):
        print "n = %d" % (n)
        pl_sub_scores = []
        lin_sub_scores = []
        for i in range(10):
            pdata = Dp.sample(n)
            qdata = Dq.sample(n)            
            PL = PluginEstimator(pdata, qdata, alpha, beta, 2)
            LIN = LinearEstimator(pdata, qdata, alpha, beta, 2);
            pl_sc = PL.eval(100)
            lin_sc = LIN.eval(100);
            print "Score(Plugin): %0.2f, Score(Linear): %0.2f" % (pl_sc, lin_sc)
            pl_sub_scores.append(pl_sc)
            lin_sub_scores.append(lin_sc)
        pl_scores.append(pl_sub_scores)
        lin_scores.append(lin_sub_scores)
    
    print scores
