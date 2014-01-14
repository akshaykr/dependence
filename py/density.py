import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

"""
Module for generating densities and related subroutines.
"""

def linear_combination(coeffs, fns, x):
    return np.sum([coeffs[i]*fns[i](x) for i in range(len(coeffs))])

class TrigDensity(object):
    """
    Multi-dimensional density object with trigonometric basis.
    The density is a linear combination of exp^{2 i pi k^Tx} for k \in \ZZ^d.
    The coefficients decay so that the function lives in the Sobolev space with smoothness s and parameter L.
    """
    def __init__(self,s,L,d):
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
        fns = [lambda x: 1]

        total = 0
        curr = [0 for i in range(self.d)]
        while total <= self.L:
            curr[np.argmin(curr)] += 1
            new_coeff = np.random.normal(0, 0.1)
            coeffs.append(new_coeff)
            fns.append(lambda x: np.exp(2j*np.pi*np.sum([curr[i]*x[i] for i in range(len(curr))])))
            total += new_coeff**2 * np.sum([curr[i]**(2*self.s) for i in range(len(curr))])
        self.coeffs = coeffs
        self.fns = fns

    def sample(self, n):
        """
        Rejection sampler for f(x) = \sum_i coeffs_i fns_i(x)
        Using uniform proposal distribution (which we hope will dominate 1/2 f(x))
        """
        x = []
        while len(x) < n:
            proposal = np.random.uniform(0, 1, [n, self.d])
            us = np.random.uniform(0,1,n)
            to_retain = [i for i in range(n) if us[i] < linear_combination(self.coeffs, self.fns, [proposal[i,0], proposal[i,1]])/2]
            if x == []:
                x = proposal[to_retain,:]
            else:
                x = np.append(x, proposal[to_retain,:], 0)
        return x[0:n,:]

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
        for i in x:
            for j in y:
                z.append(linear_combination(self.coeffs, self.fns, [i,j]))
        z = np.array(z)
        Z = np.real(z.reshape(len(x), len(y)))
        X,Y = np.meshgrid(x,y)
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
        z = []
        for i in x:
            for j in y:
                z.append(linear_combination(self.coeffs, self.fns, [i,j]))
        z = np.array(z)
        Z = np.real(z.reshape(len(x), len(y)))
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
        assert data.shape[1] == 2
        H, xedges, yedges = np.histogram2d(data[:,0], data[:,1], bins=50)
        H.shape, xedges.shape, yedges.shape = ((50, 50), (51,), (51,))
        
        extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.imshow(H, extent=extent, interpolation='nearest')
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
        fns = [lambda x: 1]
        total = 0
        curr = 2
        while total <= self.L:
            new_coeff = np.random.normal(0, 0.1)
            coeffs.append(new_coeff)
            if curr % 2 == 0:
                fns.append(lambda x: np.sqrt(2) * np.cos(np.pi*curr*x))
                total += new_coeff**2 * curr**(2*self.s)
            else:
                fns.append(lambda x: np.sqrt(2) * np.sin(np.pi*(curr-1)*x))
                total += new_coeff**2 * (curr-1)**(2*self.s)
            curr += 1
        self.coeffs = coeffs
        self.fns = fns

    def sample(self, n):
        """
        Rejection sampler for f(x) = \sum_i coeffs_i fns_i(x)
        Using uniform proposal distribution (which we hope will dominate 1/2 f(x))
        """
        x = []
        while len(x) < n:
            proposal_samples = np.random.uniform(0,1,n)
            us = np.random.uniform(0, 1, n)
            
            to_retain = [us[i] < linear_combination(self.coeffs, self.fns, proposal_samples[i])/2 for i in range(len(proposal_samples))]
            x.extend([proposal_samples[i] for i in range(len(proposal_samples)) if to_retain[i]])
        return x[0:n]

    def plot_fn(self, ax=None):
        """
        Plot the density. If supplied ax, plot on that axis. Otherwise plot and show in a new figure.
        """
        x = np.arange(0, 1, 0.01)
        y = [linear_combination(self.coeffs, self.fns, x[i]) for i in range(len(x))]
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
        if ax==None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.hist(data, bins=100)
            plt.show()
        else:
            ax.hist(data, bins=100)

if __name__=='__main__':
    fig = plt.figure()
    ax1 = fig.add_subplot(231, projection='3d')
    print "Constructing 2d Sobolev density"
    TD = TrigDensity(4,5,2)
    TD.plot_surface(ax=ax1)
    print "Rejection Sampling 2d Sobolev density"
    data = TD.sample(100000)
    ax2 = fig.add_subplot(232)
    TD.plot_fn_histogram(ax=ax2)
    ax3 = fig.add_subplot(233)
    TD.plot_data_histogram(data, ax=ax3)

    print "Constructing 1d Sobolev density"
    TD = UniTrigDensity(4,5)
    ax4 = fig.add_subplot(234)
    TD.plot_fn(ax=ax4)

    print "Rejection Sampling 1d Sobolev density"
    data = TD.sample(10000)
    ax5 = fig.add_subplot(235)
    TD.plot_data_histogram(data, ax=ax5)
    plt.show()
