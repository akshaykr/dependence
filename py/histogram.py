import numpy as np
import density, helper, kernels
import matplotlib.pyplot as plt


class Histogram(object):
    def __init__(self, data, s):
        """
        Initialize a histogram estimator with a n x d matrix of data and the appropriate smoothness s
        We need to know the smoothness to choose the correct bin width.
        TODO: only works in 1 dimension right now
        """
        self.s = s
        self.data = data
        self.n = data.shape[0]
        self.d = data.shape[1]
        self.h = 0.5*np.power(self.n*self.d, -1.0/(2*self.s+self.d))

        self.bins = np.arange(0.0, 1.3, self.h)
        self.thetas = [np.mean(np.logical_and(self.data > self.bins[i], self.data <= self.bins[i+1])) for i in range(len(self.bins)-1)]
            
    def eval(self, pts):
        m = pts.shape[0]
        out = []
        for i in range(m):
            x = pts[i,0]
            out.append(self.thetas[int(x/self.h)]/min(self.h, 1.0 - self.bins[int(x/self.h)]))
        return np.array(out)

    def kde_error(self, true_p, p_norm, fast=True):
        """
        compute the error of this estimator in ell_p^p norm. 
        """
        if fast:
            integrator = lambda x,y,z: helper.fast_integration(x,y,z)
        else:
            integrator = lambda x,y,z: helper.numeric_integration(x,y,z)
        fn_handle = lambda x: np.power(np.array(np.abs(self.eval(np.matrix(x)) - true_p.eval(np.matrix(x)).reshape(x.shape[0],)))[0,:], p_norm)
        return integrator(fn_handle, [0.1 for i in range(self.d)], [0.9 for i in range(self.d)])
