import numpy as np
import density, helper, kernels
import matplotlib.pyplot as plt


class KDE(object):

    def __init__(self, data, s):
        """
        Initialize a kernel density estimator with a n x d matrix of data and the appropriate smoothness s.
        We need to know the smoothness to choose the correct kernel and to choose the correct bandwidth.

        We correct for bias on the boundary by mirroring the data across the axis.
        """
        self.s = s
        self.m = None
        if int(self.s) == self.s:
            self.m = self.s
        else:
            self.m = int(self.s)+1
        self.data = data
        self.d = data.shape[1]
        self.n = data.shape[0]
        self.replicated_data = np.matrix(np.zeros((3**self.d*self.n, self.d)))
        for i in range(3**self.d):
            tuple = [int(i/3**j) % 3 for j in range(self.d)]
            for j in range(self.d):
                if tuple[j] == 0:
                    self.replicated_data[(self.n*i):(self.n*(i+1)), j] = -1*self.data[:,j]
                if tuple[j] == 1:
                    self.replicated_data[(self.n*i):(self.n*(i+1)), j] = self.data[:,j]
                if tuple[j] == 2:
                    self.replicated_data[(self.n*i):(self.n*(i+1)), j] = 2 - self.data[:,j]
        
        self.data = self.replicated_data
        self.h = 0.5*np.power(self.n, -1.0/(2*self.s+self.d))

        to_keep = [i for i in range(self.data.shape[0]) if np.max(np.abs(self.data[i,:] - np.matrix(0.5*np.ones((1, self.d))))) <= 0.5+self.h]
        self.data = self.replicated_data[to_keep,:]

        self.kernel = lambda x, c: kernels.kernel(x, self.m, self.h, centre=c)
        
    def eval2(self, pts):
        """
        Deprecated
        Evaluate the kernel density estimator at a set of points x.
        """
        assert False, "Deprecated"
        vals = 1.0/(self.n*self.h**self.d) * np.sum([self.kernel(self.data[j,:] - pts, None) for j in range(self.n)], axis=0).T
        return vals[0,:]

    def eval(self, pts):
        """
        Evaluate the kernel density estimator at a set of points x.
        """
        vals = 1.0/(self.n * self.h**self.d) * np.sum([self.kernel(self.data, pts[i,:]) for i in range(pts.shape[0])], axis=1)
        vals = np.maximum(vals, np.matrix(0.5*np.ones(vals.shape)))
        return vals

    def kde_error(self, true_p, p_norm, pts=1000):
        """
        compute the error of this estimator in ell_p^p norm. 
        """
        fn_handle = lambda x: np.power(np.array(np.abs(self.eval(np.matrix(x)) - true_p.eval(np.matrix(x)).reshape(x.shape[0],)))[0,:], p_norm)
        return helper.fast_integration(fn_handle, [0.0 for i in range(self.d)], [1.0 for i in range(self.d)], pts=pts)

    def kde_error2(self, true_p, p_norm):
        """
        Deprecated -- a slower way to compute the error
        """
        assert False, "Deprecated"
        coords = np.matrix(np.arange(0, 1, 0.01)).T
        vals = self.eval(coords)
        truth = true_p.eval(coords).reshape(coords.shape[0],)
        return 1.0/coords.shape[0]* np.linalg.norm(np.array(vals-truth)[0,:], ord=p_norm)**p_norm

if __name__=='__main__':
    do_one_d = True
    do_two_d = False #True
    n = 10000
    s = 2

    if do_one_d:
        print "One dimensional example"
        D = density.UniTrigDensity(s, 10)

        print "Sampling %d samples from univariate density" % (n)
        data = D.sample(n)
        
        K = KDE(data, s)
        
        pts = np.matrix(np.arange(0, 1.0, 0.001)).T
        
        print "Evaluating KDE on %d points" % (pts.shape[0]) 
        vals = K.eval(pts)

        fig = plt.figure(1, figsize=(10, 5))
        ax1 = fig.add_subplot(131)
        D.plot_fn(ax=ax1)
        plt.axis((0,1,0, 1.25));
        ax2 = fig.add_subplot(132)
        D.plot_data_histogram(data, ax=ax2)
        ax3 = fig.add_subplot(133)
        ax3.plot(pts, vals)
        plt.axis((0,1,0, 1.25));

    if do_two_d:
        print "Two dimensional example"
        D = density.TrigDensity(s, 5, 2)
        print "Sampling %d samples from 2-d density" % (n)
        data = D.sample(n)

        K = KDE(data, s)
        t = 100
        x = np.arange(0, 1, 1.0/t)
        y = np.arange(0, 1, 1.0/t)
        X,Y = np.meshgrid(x,y)
        v = np.matrix(zip(X.reshape(t**2,), Y.reshape(t**2,)))

        print "Evaluating KDE on %d points" % (len(x)*len(y)) 
        z = K.eval(v)
        z = np.array(z)
        Z = z.reshape(len(x), len(y))
        
        fig = plt.figure(2, figsize=(10,5))
        ax1 = fig.add_subplot(131)
        D.plot_fn_histogram(ax=ax1)
        ax2 = fig.add_subplot(132)
        D.plot_data_histogram(np.array(data), ax=ax2)
        ax3 = fig.add_subplot(133)
        ax3.imshow(Z, extent=[0, 1.0, 0, 1.0])

    print "Plotting"
    plt.show()
