import numpy as np
import kde, density, lattice
import itertools
from helper import *

lb = 0.1
ub = 0.9

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
        plugin_est = super(LinearEstimator, self).eval(fast=fast);

        if fast:
            integrator = lambda x, y, z: fast_integration(x,y,z)
        else:
            integrator = lambda x, y, z: numeric_integration(x,y,z)

        # C1 = - (alpha + beta) * int {p0(x)^alpha q0(x)^beta} dx
        if self.dim == 1:
          C1_fn_handle = lambda x: np.multiply( 
            np.power(self.Kp.eval(np.matrix(x)), self.alpha),
            np.power(self.Kq.eval(np.matrix(x)), self.beta) );
          l_limit = [lb];
          u_limit = [ub];

        elif self.dim == 2:
          C1_fn_handle = lambda x,y : np.multiply( 
            np.power(self.Kp.eval( np.concatenate((np.matrix(x), np.matrix(y)),
                                   1) ), self.alpha),
            np.power(self.Kq.eval( np.concatenate((np.matrix(x), np.matrix(y)),
                                   1) ), self.beta) );
          l_limit = [lb, lb];
          u_limit = [ub, ub];
        C1 = -(self.alpha + self.beta) * integrator(
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


class QuadraticEstimator(PluginEstimator):
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
        self.m = np.power(18*self.dim/s * np.power(2, 4.0*s/self.dim)* num_pdata**(-2), -float(self.dim)/(4*s+self.dim))

    def eval(self, fast=True):
        # Plugin estimator
        plugin_est = super(QuadraticEstimator, self).eval(fast=fast);

        if fast:
            integrator = lambda x, y, z: fast_integration(x,y,z)
        else:
            integrator = lambda x, y, z: numeric_integration(x,y,z)

        # C2 = 1/2(alpha(alpha-1) + alpha beta + beta(beta-1)) * int {p0(x)^alpha q0(x)^beta} dx
        if self.dim == 1:
          C2_fn_handle = lambda x: np.multiply( 
            np.power(self.Kp.eval(np.matrix(x)), self.alpha),
            np.power(self.Kq.eval(np.matrix(x)), self.beta) );
          l_limit = [lb];
          u_limit = [ub];

        elif self.dim == 2:
          C2_fn_handle = lambda x,y : np.multiply( 
            np.power(self.Kp.eval( np.concatenate((np.matrix(x), np.matrix(y)),
                                   1) ), self.alpha),
            np.power(self.Kq.eval( np.concatenate((np.matrix(x), np.matrix(y)),
                                   1) ), self.beta) );
          l_limit = [lb, lb];
          u_limit = [ub, ub];
        C2 = 0.5 * (self.alpha*(self.alpha-1) + self.alpha*self.beta + self.beta*(self.beta-1)) * integrator(C2_fn_handle, l_limit, u_limit);

        # theta^p_{2,1} = \alpha(2-\alpha+\beta/2) \EE[\phat^{\alpha-1}(X)\qhat^\beta(X)]
        theta_p_21 = self.alpha*(2-self.alpha+self.beta/2)*np.mean(
            np.multiply(np.power(self.Kp.eval(self.p_est_data), self.alpha-1),
                        np.power(self.Kq.eval(self.p_est_data), self.beta))
            )

        # theta^q_{2,1} = \beta (2 - \beta+\alpha/2) \EE[\phat^\alpha(X)\qhat^{\beta-1}(X)]
        theta_q_21 = self.beta*(2-self.beta+self.alpha/2)*np.mean(
            np.multiply(np.power(self.Kp.eval(self.q_est_data), self.alpha),
                        np.power(self.Kq.eval(self.q_est_data), self.beta-1))
            )

        # theta^{p,q}_{2,2} = 1/2 \alpha \beta \int \phat^{\alpha-1} \qhat^{\beta-1} pq
        theta_pq_22 = np.sum([
                np.mean(self.comp_exp(fn,self.p_est_data)) *
                np.mean(np.array(self.comp_exp(fn, self.q_est_data)) *
                        np.array(np.power(self.Kp.eval(self.q_est_data), self.alpha-1)) *
                        np.array(np.power(self.Kq.eval(self.q_est_data), self.beta-1)))
                for fn in lattice.lattice(self.dim, self.m)])

        theta_pq_22 *= 0.5 * self.alpha * self.beta

        # theta^p_{2,2} = 1/2 \alpha(\alpha-1) \int \phat^{\alpha-2}\qhat^\beta p^2
        ## TODO: Not using the U-statistic here because I don't think it will be that much better.
        ## Instead I'm using a related estimator as above
        theta_p_22 = 0.5*self.alpha*(self.alpha-1) * np.sum([
                np.mean(self.comp_exp(fn,self.p_est_data)) *
                np.mean(np.array(self.comp_exp(fn, self.p_est_data)) *
                        np.array(np.power(self.Kp.eval(self.p_est_data), self.alpha-2)) *
                        np.array(np.power(self.Kq.eval(self.p_est_data), self.beta)))
                for fn in lattice.lattice(self.dim, self.m)])
        

        # theta^q_{2,2} = 1/2 \beta(\beta-1) \int \phat^{\alpha}\qhat^{\beta-2} q^2
        theta_q_22 = 0.5 * self.beta*(self.beta-1) * np.sum([
                np.mean(self.comp_exp(fn,self.q_est_data)) *
                np.mean(np.array(self.comp_exp(fn, self.q_est_data)) * 
                        np.array(np.power(self.Kp.eval(self.q_est_data), self.alpha)) * 
                        np.array(np.power(self.Kq.eval(self.q_est_data), self.beta-2)))
                for fn in lattice.lattice(self.dim, self.m)])

        return np.real(plugin_est + theta_p_21 + theta_q_21 + theta_p_22 + theta_q_22 + theta_pq_22 + C2);

    def quad_term_estimator(self, fn, data1, data2):
        np.mean(self.comp_exp(k, data1) * np.mean(np.array))

    def comp_exp(self, fn, x):
        return np.exp(2j*np.pi*np.matrix(fn)*np.matrix(x).T)

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
