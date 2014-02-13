import numpy as np
import kde, density, lattice
# TODO. WTF IS THIS * Import
from helper import (
    fast_integration,
    numeric_integration,
)

lb = 0.0
ub = 1.0

class PluginEstimator(object):
    """
    Plugin estimator for functionals of the form \int
    p^\alpha(x)q^\beta(x). We construct kernel density estimators for
    p,q$ and then perform numerical integration to evalute \int \phat^\alpha(x) \qhat^\beta(x)
    """
    def __init__(self, pdata, qdata, alpha, beta, s):
        """
        Constructor. Takes in
        pdata -- data drawn from p
        qdata -- data drawn from q
        alpha -- exponent for p
        beta -- exponent for q
        s -- smoothness (which we require) of the two densities.
        """
        self.Kp = kde.KDE(pdata, s)
        self.Kq = kde.KDE(qdata, s)
        self.dim = pdata.shape[1]
        self.alpha = alpha
        self.beta = beta
        self.s = s

    def eval(self, fast=True):
        """
        Evaluate the estimator by performing numeric integration on the n^d grid of points.
        Currently we ignore the fast parameter -- we are always doing fast numeric integration.
        """
        val = fast_integration(lambda x: np.multiply(
                np.power(self.Kp.eval(np.matrix(x)), self.alpha),
                np.power(self.Kq.eval(np.matrix(x)), self.beta)),
                               [lb for i in range(self.dim)], [ub for i in range(self.dim)])
        return val


class LinearEstimator(PluginEstimator):
    """
    Linear estimator for functionals of the form \int p^\alpha(x)q^\beta(x).
    We construct kernel density estimators for p,q on the first half of the sample and then evaluate three terms:
    1. (1 - \alpha - beta) \int \phat^\alpha \qhat^\beta.
    2. alpha \int \phat^{\alpha-1} \qhat^\beta p.
    3. beta \int \phat^\alpha \qhat^{\beta-1} q.
    The last two terms are evaluated via sample means on the second half of the sample.
    We perform numerical integration to evalute \int \phat^\alpha(x) \qhat^\beta(x)
    """
    def __init__(self, pdata, qdata, alpha, beta, s):
        """
        Constructor. Takes in
        pdata -- data drawn from p
        qdata -- data drawn from q
        alpha -- exponent for p
        beta -- exponent for q
        s -- smoothness (which we require) of the two densities.
        """
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

    def eval(self, fast=True):
        """
        Evaluate the estimator.
        Currently we ignore the fast parameter -- we are always doing fast numeric integration.
        """
        # Plugin estimator
        c1 = 1 - self.alpha - self.beta
        if c1 != 0:
            print "Calling plugin estimator"
            plugin_est = c1 * super(LinearEstimator, self).eval(fast=fast);
        else:
            plugin_est = 0.0

        # theta^p_{1,1} = alpha * E_{x~p}[p0(x)^(alpha-1) * q0(x)^beta]
        theta_p_11 = self.alpha * self.linear_term(lambda x:
                                                       np.multiply(np.power(self.Kp.eval(x), self.alpha-1),
                                                                   np.power(self.Kq.eval(x), self.beta)),
                                                   self.p_est_data)

        # theta^q_{1,1} = beta * E_{x~q}[p0(x)^alpha * q0(x)^{beta-1}]
        theta_q_11 = self.beta * self.linear_term(lambda x:
                                                       np.multiply(np.power(self.Kp.eval(x), self.alpha),
                                                                   np.power(self.Kq.eval(x), self.beta-1)),
                                                   self.q_est_data)
        # Return plugin + theta_p_11 + theta_q_11
        return plugin_est + theta_p_11 + theta_q_11;

    def linear_term(self, fn, data):
        """
        Helper routine for computing the mean of a function evaluated on a data set.
        """
        return np.mean(fn(data))

class QuadraticEstimator(PluginEstimator):
    """
    Quadratic estimator for functionals of the form \int p^\alpha(x)q^\beta(x).
    We construct kernel density estimators for p,q on the first half of the sample and then evaluate a bunch of terms.
    See the paper for details
    We perform numerical integration to evalute several terms.
    """
    def __init__(self, pdata, qdata, alpha, beta, s):
        """
        Constructor. Takes in
        pdata -- data drawn from p
        qdata -- data drawn from q
        alpha -- exponent for p
        beta -- exponent for q
        s -- smoothness (which we require) of the two densities.
        """
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
        self.m = 2*np.power(18*self.dim/s * np.power(2, 4.0*s/self.dim)* num_pdata**(-2), -float(self.dim)/(4*s+self.dim))

    def eval(self, fast=True):
        """
        Evaluate the estimator.
        Currently we ignore the fast parameter -- we are always doing fast numeric integration.
        """
        # Plugin estimator
        c2 = 1 - 1.5*(self.alpha+self.beta) + 0.5*(self.alpha+self.beta)**2
        if c2 != 0:
            plugin_est = c2 * super(QuadraticEstimator, self).eval(fast=fast);
        else:
            plugin_est = 0.0

        ## estimate the linear terms
        # theta^p_{2,1} = \alpha(2-\alpha-\beta) \EE[\phat^{\alpha-1}(X)\qhat^\beta(X)]
        theta_p_21 = self.alpha*(2-self.alpha-self.beta)*np.mean(
            np.multiply(np.power(self.Kp.eval(self.p_est_data), self.alpha-1),
                        np.power(self.Kq.eval(self.p_est_data), self.beta))
            )

        # theta^q_{2,1} = \beta (2 -\beta-\alpha) \EE[\phat^\alpha(X)\qhat^{\beta-1}(X)]
        theta_q_21 = self.beta*(2-self.beta-self.alpha)*np.mean(
            np.multiply(np.power(self.Kp.eval(self.q_est_data), self.alpha),
                        np.power(self.Kq.eval(self.q_est_data), self.beta-1))
            )

        ## Estimate the quadratic terms
        # theta^{p,q}_{2,2} = \alpha \beta \int \phat^{\alpha-1} \qhat^{\beta-1} pq
        theta_pq_22 = self.bilinear_term_fast(lambda x: 
                                              np.array(np.power(self.Kp.eval(x), self.alpha-1)) *
                                              np.array(np.power(self.Kq.eval(x), self.beta-1)),
                                              self.p_est_data,
                                              self.q_est_data)

        theta_pq_22 = self.alpha * self.beta * theta_pq_22

        # theta^p_{2,2} = 1/2 \alpha(\alpha-1) \int \phat^{\alpha-2}\qhat^\beta p^2
        theta_p_22 = 0.5*self.alpha*(self.alpha-1) * self.quad_term_fast(
            lambda x: np.array(np.power(self.Kp.eval(x), self.alpha-2)) * np.array(np.power(self.Kq.eval(x), self.beta)),
            self.p_est_data)

        # theta^q_{2,2} = 1/2 \beta(\beta-1) \int \phat^{\alpha}\qhat^{\beta-2} q^2
        theta_q_22 = 0.5 * self.beta*(self.beta-1) * self.quad_term_fast(
            lambda x: np.array(np.power(self.Kp.eval(x), self.alpha)) * np.array(np.power(self.Kq.eval(x), self.beta-2)),
            self.q_est_data)

        return np.real(plugin_est + theta_p_21 + theta_q_21 + theta_p_22 + theta_q_22 + theta_pq_22);

    def bilinear_term_slow(self, fn, data1, data2):
        """
        Deprected
        A slow routine for estimating the bilinear term. This iterates
        over all of the basis elements and then over all of the points
        in the sample.
        fn -- a lambda expression for the function \psi.
        data1 -- data from the distribution p
        data2 -- data from the distribution q

        Returns the value of the estimator
        """
        assert False, "Deprected"
        n1 = data1.shape[0]
        n2 = data2.shape[0]
        total = 0.0
        for k in lattice.lattice(self.dim, self.m):
            for i in range(n1):
                for j in range(n2):
                    total += self.comp_exp(k, data1[i,:])*self.comp_exp(k, data2[j,:])*fn(data2[j,:])
        return np.real(total[0,0]/(n1*n2))

    def bilinear_term_fast(self, fn, data1, data2):
        """
        Faster routine for estimating the bilinear term.
        fn -- a lambda expression for the function \psi.
        data1 -- data from the distribution p
        data2 -- data from the distribution q

        Returns the value of the estimator

        Note: This shouldn't be called externally.
        """
        total = np.sum([
                np.mean(self.comp_exp(k,data1)) *
                np.mean(np.array(self.comp_exp(k, data2)) *
                        np.array(fn(data2)))
                for k in lattice.lattice(self.dim, self.m)])
        return np.real(total)

    def quad_term_slow(self, fn, data):
        """
        Deprected
        Slow routine for estimating the quadratic terms
        fn -- a lambda expression for the function \psi.
        data1 -- data from the distribution we are trying to estimate

        Returns the value of the estimator
        """
        assert False, "Deprecated"
        n = data.shape[0]
        total = 0.0
        for k in lattice.lattice(self.dim, self.m):
            for i in range(n):
                for j in range(n):
                    if j != i:
                        total += self.comp_exp(k, data[i,:])*self.comp_exp(k, data[j,:])*fn(data[j,:])
        term2 = 0.0
        for k in lattice.lattice(self.dim, self.m):
            for kp in lattice.lattice(self.dim, self.m):
                bi = fast_integration(lambda x: np.array(self.comp_exp(k,x))*np.array(self.comp_exp(kp,x))*np.array(fn(x)),
                                      [0 for t in range(self.dim)], [1 for t in range(self.dim)])
                for i in range(n):
                    for j in range(n):
                        if j != i:
                            term2 += bi*self.comp_exp(k, data[i,:])*self.comp_exp(kp, data[j,:])

        return np.real(2.0*total/(n*(n-1)) - 1.0*term2/(n*(n-1)))

    def quad_term_fast(self, fn, data):
        """
        Fast routine for estimating the quadratic term
        fn -- a lambda expression for the function \psi.
        data1 -- data from the distribution we are trying to estimate

        Returns the value of the estimator

        Note: this shouldn't be called externally
        """
        n = data.shape[0]
        total = 0.0
        for k in lattice.lattice(self.dim, self.m):
            sub1 = np.array(self.comp_exp(k, data))
            sub2 = sub1*np.array(fn(data))
            total += np.sum((np.sum(sub2) - sub2) * sub1)

        term2 = 0.0
        for k in lattice.lattice(self.dim, self.m):
            for kp in lattice.lattice(self.dim, self.m):
                bi = fast_integration(lambda x: np.array(self.comp_exp(k,x))*np.array(self.comp_exp(kp,x))*np.array(fn(x)),
                                      [0.0 for t in range(self.dim)], [1.0 for t in range(self.dim)], pts=100)
                sub1 = np.array(self.comp_exp(k, data))
                sub2 = np.array(self.comp_exp(kp, data))
                term2 += np.sum(bi* sub1* (np.sum(sub2) - sub2))
#         print (np.real(2.0*total/(n*(n-1))), np.real(term2/(n*(n-1))))
        return np.real(2.0*total/(n*(n-1))) - np.real(term2/(n*(n-1)))

    def comp_exp(self, fn, x):
        """
        evaluate the complex exponential with parameter fn (1-d vector) on the set of points x (n-by-d matrix)
        """
        return np.exp(2j*np.pi*np.matrix(fn)*np.matrix(x).T).T

class L2Estimator(object):
    """
    Estimator for the L_2^2 divergence.
    """
    def __init__(self, pdata, qdata, s):
        """
        Constructor.
        pdata -- data from p
        qdata -- data from q
        s -- smoothness of the densities. We need to know the smoothness.
        """
        self.pdata = pdata
        self.qdata = qdata
        self.dim = pdata.shape[1]
        self.n = pdata.shape[0]
        self.s = s
        self.m = 2*np.power(18*self.dim/self.s * np.power(2, 4.0*self.s/self.dim)* self.n**(-2), -float(self.dim)/(4*self.s+self.dim))

    def eval(self, fast=True):
        """
        Evaluate the estimator on the data.
        See the paper for details about the estimator.
        """
        T1 = 0.0
        T2 = 0.0
        ## T1 = \int p^2. T2 = \int q^2, T3 = -2 \int pq
        for k in lattice.lattice(self.dim, self.m):
            T1 += 1.0/self.pdata.shape[0]**2 *(np.sum(np.array(self.comp_exp(k, self.pdata)))**2 - np.sum(np.array(self.comp_exp(k,self.pdata))**2))
            T2 += 1.0/self.qdata.shape[0]**2 *(np.sum(np.array(self.comp_exp(k, self.qdata)))**2 - np.sum(np.array(self.comp_exp(k,self.qdata))**2))

        T3 = -2.0 * np.sum([
                np.mean(self.comp_exp(k,self.pdata)) *
                np.mean(np.array(self.comp_exp(k, self.qdata)))
                for k in lattice.lattice(self.dim, self.m)])
        return np.real(T1 + T2 + T3)

    def comp_exp(self, fn, x):
        """
        evaluate the complex exponential with parameter fn (1-d vector) on the set of points x (n-by-d matrix)
        """
        return np.exp(2j*np.pi*np.matrix(fn)*np.matrix(x).T).T

class Truth(object):
    def __init__(self, p, q, alpha, beta):
        self.p = p
        self.q = q
        self.alpha = alpha
        self.beta = beta
        self.dim = self.p.d

    def eval(self, fast=True, pts=1000):
        """
        Evaluate the true divergence by performing numeric integration on the n^d grid of points.
        """
        if not fast:
            print "using slow integration"
            if self.dim == 1:
                val = numeric_integration(lambda x: np.multiply(
                        np.power(self.p.eval(np.matrix(x)), self.alpha),
                        np.power(self.q.eval(np.matrix(x)), self.beta)),
                                          [lb], [ub])
            if self.dim == 2:
                val = numeric_integration(lambda x,y: np.multiply(
#                         np.power(self.p.eval(np.concatenate(np.matrix(x), np.matrix(y))), self.alpha),
#                         np.power(self.q.eval(np.concatenate(np.matrix(x), np.matrix(y))), self.beta)),
                        np.power(self.p.eval(np.concatenate(np.matrix(x), np.matrix(y))), self.alpha),
                        np.power(self.q.eval(np.concatenate(np.matrix(x), np.matrix(y))), self.beta)),
                                          [lb, lb], [ub, ub])

        if fast:
            print "using fast integration"
            if self.dim == 1:
                val = fast_integration(lambda x: np.multiply(
                        np.power(self.p.eval(np.matrix(x)), self.alpha),
                        np.power(self.q.eval(np.matrix(x)), self.beta)),
                                          [lb], [ub], pts=pts)
            if self.dim == 2:
                val = fast_integration(lambda x: np.multiply(
                        np.power(self.p.eval(np.matrix(x)), self.alpha),
                        np.power(self.q.eval(np.matrix(x)), self.beta)),
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
