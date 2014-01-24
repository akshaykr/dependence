import numpy as np
import kde, density, lattice
import itertools
from helper import *

lb = 0.0
ub = 1.0

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
                                       [lb], [ub])
            else:
                val = numeric_integration(lambda x: np.multiply(
                        np.power(self.Kp.eval(np.matrix(x)), self.alpha), 
                        np.power(self.Kq.eval(np.matrix(x)), self.beta)),
                                          [lb], [ub])
        if self.dim == 2:
            if fast:
                val = fast_integration(lambda x: np.multiply(
                        np.power(self.Kp.eval(np.matrix(x)), self.alpha),
                        np.power(self.Kq.eval(np.matrix(x)), self.beta)),
                                          [lb,lb], [ub,ub])
            else:
                val = numeric_integration(lambda x,y: np.multiply(
                        np.power(self.Kp.eval(np.concatenate((np.matrix(x), np.matrix(y)),
                                                             1)), self.alpha),
                        np.power(self.Kq.eval(np.concatenate((np.matrix(x), np.matrix(y)),
                                                             1)), self.beta)),
                                          [lb,lb], [ub,ub])
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

    def eval(self, fast=False):
        # Plugin estimator
        c1 = 1 - self.alpha - self.beta
        if c1 != 0:
            plugin_est = c1 * super(LinearEstimator, self).eval(fast=fast);
        else:
            plugin_est = 0.0

        if fast:
            integrator = lambda x, y, z: fast_integration(x,y,z)
        else:
            integrator = lambda x, y, z: numeric_integration(x,y,z)

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
        return plugin_est + theta_p_11 + theta_q_11;


class QuadraticEstimator(PluginEstimator):
    def __init__(self, pdata, qdata, alpha, beta, s):
        # Shuffle the data and split it into 2 for density estimation and
        # estimating the other parts.
#         np.random.shuffle(pdata);
#         np.random.shuffle(qdata);
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
        self.m = 2*np.power(18*self.dim/s * np.power(2, 4.0*s/self.dim)* num_pdata**(-2), -float(self.dim)/(4*s+self.dim))

    def eval(self, fast=True):
        # Plugin estimator
        c2 = 1 - 1.5*(self.alpha+self.beta) + 0.5*(self.alpha+self.beta)**2
        if c2 != 0:
            plugin_est = c2 * super(QuadraticEstimator, self).eval(fast=fast);
        else:
            plugin_est = 0.0

        if fast:
            integrator = lambda x, y, z: fast_integration(x,y,z)
        else:
            integrator = lambda x, y, z: numeric_integration(x,y,z)

        # theta^p_{2,1} = \alpha(2-\alpha-\beta) \EE[\phat^{\alpha-1}(X)\qhat^\beta(X)]
        theta_p_21 = self.alpha*(2-self.alpha-self.beta)*np.mean(
            np.multiply(np.power(self.Kp.eval(self.p_est_data), self.alpha-1),
                        np.power(self.Kq.eval(self.p_est_data), self.beta))
            )

        # theta^q_{2,1} = \beta (2 -\beta-\alpha) \EE[\phat^\alpha(X)\qhat^{\beta-1}(X)]
        theta_q_21 = self.beta*(2-self.beta-self.alpha/2)*np.mean(
            np.multiply(np.power(self.Kp.eval(self.q_est_data), self.alpha),
                        np.power(self.Kq.eval(self.q_est_data), self.beta-1))
            )

        # theta^{p,q}_{2,2} = \alpha \beta \int \phat^{\alpha-1} \qhat^{\beta-1} pq
        theta_pq_22 = self.bilinear_term_fast(lambda x: 
                                              np.array(np.power(self.Kp.eval(x), self.alpha-1)) *
                                              np.array(np.power(self.Kq.eval(x), self.beta-1)),
                                              self.p_est_data,
                                              self.q_est_data)

        theta_pq_22 = self.alpha * self.beta * theta_pq_22

        # theta^p_{2,2} = 1/2 \alpha(\alpha-1) \int \phat^{\alpha-2}\qhat^\beta p^2
        theta_p_22 = 0.5*self.alpha*(self.alpha-1) * self.quad_term_slow(
            lambda x: np.array(np.power(self.Kp.eval(x), self.alpha-2)) * np.array(np.power(self.Kq.eval(x), self.beta)),
            self.p_est_data)

        # theta^q_{2,2} = 1/2 \beta(\beta-1) \int \phat^{\alpha}\qhat^{\beta-2} q^2
        theta_q_22 = 0.5 * self.beta*(self.beta-1) * self.quad_term_slow(
            lambda x: np.array(np.power(self.Kp.eval(x), self.alpha)) * np.array(np.power(self.Kq.eval(x), self.beta-2)),
            self.q_est_data)

        return np.real(plugin_est + theta_p_21 + theta_q_21 + theta_p_22 + theta_q_22 + theta_pq_22);

    def bilinear_term_slow(self, fn, data1, data2):
         n1 = data1.shape[0]
         n2 = data2.shape[0]
         total = 0.0
         for k in lattice.lattice(self.dim, self.m):
             for i in range(n1):
                 for j in range(n2):
                     total += self.comp_exp(k, data1[i,:])*self.comp_exp(k, data2[j,:])*fn(data2[j,:])
         return np.real(total[0,0]/(n1*n2))

    def bilinear_term_fast(self, fn, data1, data2):
        total = np.sum([
                np.mean(self.comp_exp(k,data1)) *
                np.mean(np.array(self.comp_exp(k, data2)) *
                        np.array(fn(data2)))
                for k in lattice.lattice(self.dim, self.m)])
        return np.real(total)

    def quad_term_slow(self, fn, data):
        n = data.shape[0]
        total = 0.0
        for k in lattice.lattice(self.dim, self.m):
            for i in range(n):
                for j in range(n):
                    if j != i:
                        total += self.comp_exp(k, data[i,:])*self.comp_exp(k, data[j,:])*fn(data[j,:])
        next = 0.0
        for k in lattice.lattice(self.dim, self.m):
            for kp in lattice.lattice(self.dim, self.m):
                bi = fast_integration(lambda x: np.array(self.comp_exp(k,x))*np.array(self.comp_exp(kp,x))*np.array(fn(x)), 
                                      [0 for t in range(self.dim)], [1 for t in range(self.dim)])
                for i in range(n):
                    for j in range(n):
                        if j != i:
                            next += bi*self.comp_exp(k, data[i,:])*self.comp_exp(kp, data[j,:])

        print (np.real(2.0*total/(n*(n-1))), np.real(next/(n*(n-1))))
        return np.real(2.0*total/(n*(n-1)) - 1.0*next/(n*(n-1)))

    def comp_exp(self, fn, x):
        return np.exp(2j*np.pi*np.matrix(fn)*np.matrix(x).T).T

class Truth(object):
    def __init__(self, p, q, alpha, beta):
        self.p = p
        self.q = q
        self.alpha = alpha
        self.beta = beta
        self.dim = self.p.d

    def eval(self, fast=False):
        """
        Evaluate the true divergence by performing numeric integration on the n^d grid of points.
        """
        if fast:
            if self.dim == 1:
                val = fast_integration(lambda x: np.multiply(
                        np.power(self.p.eval(np.matrix(x)), self.alpha), 
                        np.power(self.q.eval(np.matrix(x)), self.beta)),
                                          [lb], [ub])
            if self.dim == 2:
                val = fast_integration(lambda x: np.multiply(
                        np.power(self.p.eval(np.matrix(x)), self.alpha),
                        np.power(self.q.eval(np.matrix(x)), self.beta)),
                                          [lb,lb], [ub,ub])
        if not fast:
            if self.dim == 1:
                val = numeric_integration(lambda x: np.multiply(
                        np.power(self.p.eval(np.matrix(x)), self.alpha), 
                        np.power(self.q.eval(np.matrix(x)), self.beta)),
                                          [lb], [ub])
            if self.dim == 2:
                val = numeric_integration(lambda x,y: np.multiply(
                        np.power(self.p.eval(np.concatenate((np.matrix(x), np.matrix(y)),
                                                            1)), self.alpha),
                        np.power(self.q.eval(np.concatenate((np.matrix(x), np.matrix(y)),
                                                            1)), self.beta)),
                                          [lb,lb], [ub,ub])
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
    
    print (pl_scores, lin_scores)
