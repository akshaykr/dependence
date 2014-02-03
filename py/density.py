import lattice
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import cm

"""
Module for generating sobolev densities and related subroutines.
"""

def linear_combination(coeffs, fns, x):
    """
    Computes the linear combination of complex exponential \sum_{t=1}^T coeffs[t] e^{2\pi i fns[t]^Tx}.
    coeffs should be a list of floats
    fns should be a list of d-dimensional row vectors
    x should be a n-by-d matrix of points to evaluate the function on.

    returns: an n-by-1 matrix of f(x).
    """
    z = np.sum([coeffs[i]*np.exp(2j*np.pi*fns[i]*np.matrix(x).T)[0,:] for i in range(len(coeffs))], axis=0)
    return np.abs(z).T

class TrigDensity(object):
    """
    Multi-dimensional density object with trigonometric basis.
    The density is a linear combination of exp^{2 i pi k^Tx} for k \in \ZZ^d.
    The coefficients decay so that the function lives in the Sobolev space with smoothness s and parameter L.
    """
    def __init__(self, s, L, d):
        """
        Init, s = sobolev smoothness, L = sobolev length, d = dimension
        """
        self.s = s
        self.L = L
        self.d = d
        self.generate_density()

    def generate_density(self):
        """
        Construct a d-dimensional density in the Sobolev class W(s,L)
        under the trigonometric basis \phi_k(x) = e^{2i \pi k^Tx}
        """
        coeffs = [1]
        fns = [np.matrix([0 for i in range(self.d)])]

        total = 0
        f = lattice.lattice(self.d, limit=None)
        f.next()
        while total <= self.L :
            curr = np.matrix(f.next())
            new_coeff = np.random.uniform(-0.1, 0.1)
            coeffs.append(new_coeff)
            fns.append(curr)
            total += new_coeff**2 * np.sum([np.abs(curr[i,0])**(2*self.s) for i in range(curr.shape[0])])
        self.coeffs = coeffs[0:len(coeffs)-1]
        self.fns = fns[0:len(fns)-1]

    def eval(self, pts):
        """
        Evaluate the density at the pts passed in. Pts should be a n-by-d matrix.
        returns a n-by-1 matrix of the evaluations.
        """
        vals =  linear_combination(self.coeffs, self.fns, pts)
        return np.matrix(vals)

    def sample(self, n):
        """
        Rejection sampler for f(x) = \sum_i coeffs_i fns_i(x)
        Using uniform proposal distribution (which we hope will dominate 1/2 f(x))
        input: n -- number of points to draw
        output: n-by-d matrix of samples drawn i.i.d. from the density.
        """
        x = []
        while len(x) < n:
            proposal = np.matrix(np.random.uniform(0, 1, [n, self.d]))
            us = np.matrix(np.random.uniform(0,1,n)).T
            to_retain = np.array(us < self.eval(proposal)/2.0).reshape(n,)
            if x == []:
                x = proposal[to_retain,:]
            else:
                x = np.append(x, proposal[to_retain,:], 0)
        return np.matrix(x[0:n,:])

    def plot_surface(self, ax=None):
        """
        For a 2-d density, generate a 2-d plot of the density contours.
        If supplied ax, plot on that axis. Otherwise plot and show in a new figure.
        If ax is supplied, it should be an Axis3D object.
        """
        assert self.d == 2, "plot_surface is only available for 2-d densities"
        z = []
        x = np.arange(0, 1, 0.01)
        y = np.arange(0, 1, 0.01)
        X,Y = np.meshgrid(x,y)
        v = np.matrix(zip(X.reshape(len(x)*len(y),), Y.reshape(len(x)*len(y),)))
        z = np.array(self.eval(v))
        Z = np.abs(z.reshape(len(x), len(y)))
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.plot_surface(X,Y,Z, cmap=cm.cool)
            plt.show()
        else:
            ax.plot_surface(X,Y,Z, cmap=cm.cool)
        return (X,Y,Z)

    def plot_fn_histogram(self, ax=None):
        """
        For a 2-d density, plot a top-down histogram of the density contours. This is just a colored matrix.
        If supplied ax, plot on that axis. Otherwise plot and show in a new figure.
        """
        assert self.d == 2, "plot_fn_histogram is only available for 2-d densities"
        x = np.arange(0, 1, 0.01)
        y = np.arange(0, 1, 0.01)
        X,Y = np.meshgrid(x,y)
        v = np.matrix(zip(X.reshape(len(x)*len(y),), Y.reshape(len(x)*len(y),)))
        z = np.array(self.eval(v))
        Z = np.abs(z.reshape(len(x), len(y)))
        if ax == None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.imshow(Z, extent=[0, 1.0, 0, 1.0], interpolation='nearest')
            plt.show()
        else:
            ax.imshow(Z, extent=[0, 1.0, 0, 1.0], interpolation='nearest')

    def plot_data_histogram(self, data, ax=None):
        """
        For samples from a 2-d density, plot the top-down histogram of the empirical density contours. This is a colored matrix.
        If supplied ax, plot on that axis. Otherwise plot and show in a new figure.
        """
        data = np.array(data)
        assert data.shape[1] == 2
        H, xedges, yedges = np.histogram2d(data[:,0], data[:,1], bins=50, normed=True)
        H.shape, xedges.shape, yedges.shape = ((50, 50), (51,), (51,))

        extent = [yedges[-1], yedges[0], xedges[0], xedges[-1]]
        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.imshow(H.T, extent=extent, interpolation='nearest')
            plt.show()
        else:
            ax.imshow(H, extent=extent, interpolation='nearest')

class UniTrigDensity(object):
    """
    Univariate trigonometric density object.
    The density belongs to the Sobolev space W(s,L)
    """
    def __init__(self,s,L,d=1):
        """
        Construct a univariate density. s = Sobolev smoothness parametr, L, Sobolev length.
        """
        assert d==1, "Univariate density must have d=1"
        self.s = s
        self.L = L
        self.d = d
        self.generate_density()

    def generate_density(self):
        """
        Construct a trigonometric density in one dimension
        """
        coeffs = [1]
        fns = [np.matrix([[0]])]
        total = 0
        curr = 2
        while total <= self.L:
            new_coeff = np.random.uniform(-0.1, 0.1)
            coeffs.append(new_coeff)
            fns.append(np.matrix([[curr]]))
            total += new_coeff**2 * np.abs(curr)**(2*self.s)
#             if curr % 2 == 0:
#                 fns.append(lambda x: np.sqrt(2) * np.cos(np.pi*curr*x))

#             else:
#                 fns.append(lambda x: np.sqrt(2) * np.sin(np.pi*(curr-1)*x))
#                 total += new_coeff**2 * (curr-1)**(2*self.s)
            curr += 1
        self.coeffs = coeffs[0:len(coeffs)-1]
        self.fns = fns[0:len(fns)-1]

    def eval(self, pts):
        """
        Evaluate the density at the pts passed in. Pts should be a n-by-d matrix.
        returns a n-by-1 matrix of the evaluations.
        """
        vals =  linear_combination(self.coeffs, self.fns, pts)
        return np.matrix(vals)

    def sample(self, n):
        """
        Rejection sampler for f(x) = \sum_i coeffs_i fns_i(x)
        Using uniform proposal distribution (which we hope will dominate 1/2 f(x))
        """
        x = []
        while len(x) < n:
            proposal = np.matrix(np.random.uniform(0, 1, [n, 1]))
            us = np.matrix(np.random.uniform(0,1,n)).T
            to_retain = np.array(us < self.eval(proposal)/2.0).reshape(n,)
            if x == []:
                x = proposal[to_retain,:]
            else:
                x = np.append(x, proposal[to_retain,:], 0)
        return np.matrix(x[0:n])

    def sample2(self, n, bins=1000):
        """
        deprecated
        """
        assert False, "Call to a deprecated routine"
        x = np.matrix(np.arange(0, 1, 1.0/bins)).T
        vals = self.eval(x)
        dist = np.array(vals).T/np.sum(vals)
        draws = np.random.multinomial(1, dist[0,:], n)
        uni = np.random.uniform(0, 1.0/bins, (n,1))
        return np.matrix(draws)*x + uni

    def plot_fn(self, ax=None):
        """
        Plot the density. If supplied ax, plot on that axis. Otherwise plot and show in a new figure.
        """
        x = np.matrix(np.arange(0, 1, 0.01)).T
        y = self.eval(x)
        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.plot(x,y)
            plt.show()
        else:
            ax.plot(x,y)
        return (x,y)

    def plot_data_histogram(self, data, ax=None):
        """
        Plot histogram of data.
        If supplied ax, plot on that axis. Otherwise plot and show in a new figure.
        """
        data = np.array(data)
        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(data, bins=100, normed=True)
            plt.show()
        else:
            ax.hist(data, bins=100)

class QuadraticDensity(object):
    """
    This is the distribution f(x) = 5/4 - (1/2-x)^2
    """
    def __init__(self):
        self.d = 1

    def sample(self, n):
        """
        Rejection sampler for f(x) = \sum_i coeffs_i fns_i(x)
        Using uniform proposal distribution (which we hope will dominate 1/2 f(x))
        """
        x = []
        while len(x) < n:
            proposal = np.matrix(np.random.uniform(0, 1, [n, 1]))
            us = np.matrix(np.random.uniform(0,1,n)).T
            to_retain = np.array(us < self.eval(proposal)/2.0).reshape(n,)
            if x == []:
                x = proposal[to_retain,:]
            else:
                x = np.append(x, proposal[to_retain,:], 0)
        return np.matrix(x[0:n])

    def eval(self, pts):
        return 5.0/4 - np.power(0.5-pts, 2)


if __name__=='__main__':
    fig = plt.figure()
    ax1 = fig.add_subplot(231, projection='3d')
    print "Constructing 2d Sobolev density"
    TD = TrigDensity(8,2,2)
    TD.plot_surface(ax=ax1)
    print "Rejection Sampling 2d Sobolev density"
    data = TD.sample(50000)
    ax2 = fig.add_subplot(232)
    TD.plot_fn_histogram(ax=ax2)
    ax3 = fig.add_subplot(233)
    TD.plot_data_histogram(data, ax=ax3)

    print "Constructing 1d Sobolev density"
    TD = UniTrigDensity(8,5)
    ax4 = fig.add_subplot(234)
    TD.plot_fn(ax=ax4)

    print "Rejection Sampling 1d Sobolev density"
    data = TD.sample(10000)
    ax5 = fig.add_subplot(235)
    TD.plot_data_histogram(data, ax=ax5)
    plt.show()
