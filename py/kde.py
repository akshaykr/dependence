import numpy as np
import density
import kernels
import matplotlib.pyplot as plt

class KDE(object):

    def __init__(self, data, s):
        """
        Initialize a kernel density estimator with a n x d matrix of data and the appropriate smoothness s.
        We need to know the smoothness to choose the correct kernel and to choose the correct bandwidth.
        """
        self.s = s
        self.m = None
        if int(self.s) == self.s:
            self.m = self.s-1
        else:
            self.m = int(self.s)
        self.data = data
        self.n = data.shape[0]
        self.d = data.shape[1]
        self.h = self.n**(-1.0/(2*self.s+self.d))

        ## TODO: what order am I supposed to use? s or \lfloor s \rfloor
        self.kernel = kernels.kernel
        
    def eval(self, pts):
        """
        Evaluate the kernel density estimator at a set of points x.
        """
        t = pts.shape[0]
        vals = []
        for i in range(t):
            vals.append(1.0/self.h**self.d*np.mean([self.kernel((self.data[j,:]-pts[i,:])/self.h, self.m, self.h) for j in range(self.n)]))
        return vals
    

if __name__=='__main__':
    do_one_d = True
    do_two_d = False
    n = 5000
    s = 2

    if do_one_d:
        print "One dimensional example"
        D = density.UniTrigDensity(s, 10)

        print "Sampling %d samples from univariate density" % (n)
        data = D.sample(n)
        
        K = KDE(data, s)
        
        pts = np.matrix(np.arange(0, 1.0, 0.01)).T
        
        print "Evaluating KDE on %d points" % (pts.shape[0]) 
        vals = K.eval(pts)
        
        fig = plt.figure(1, figsize=(10, 5))
        ax1 = fig.add_subplot(131)
        D.plot_fn(ax=ax1)
        ax2 = fig.add_subplot(132)
        D.plot_data_histogram(data, ax=ax2)
        ax3 = fig.add_subplot(133)
        ax3.plot(pts, vals)

    if do_two_d:
        print "Two dimensional example"
        D = density.TrigDensity(s, 5, 2)
        print "Sampling %d samples from 2-d density" % (n)
        data = D.sample(n)

        K = KDE(data, s)
        x = np.arange(0, 1, 0.1)
        y = np.arange(0, 1, 0.1)
        z = []
        
        print "Evaluating KDE on %d points" % (len(x)*len(y)) 
        for i in x:
            for j in y:
                z.append(K.eval(np.matrix([[i,j]])))
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
