import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm

"""
Module for generating densities and related subroutines.
"""

def trig_density(s, L, d):
    """
    Construct a d-dimensional density in the Sobolev class W(s,L)
    under the trigonometric basis \phi_k(x) = e^{2i \pi k^Tx}
    """
    coeffs = [1]
    fns = [lambda x: 1]


    total = 0
    curr = [0 for i in range(d)]
    while total <= L:
        curr[np.argmin(curr)] += 1
        new_coeff = np.random.normal(0, 0.1)
        coeffs.append(new_coeff)
        fns.append(lambda x: np.exp(2j*np.pi*np.sum([curr[i]*x[i] for i in range(len(curr))])))
        total += new_coeff**2 * np.sum([curr[i]**(2*s) for i in range(len(curr))])
        
    return [coeffs, fns]

def uni_trig_density(beta, L = 10):
    """
    Construct a trigonometric density in one dimension with
    sobolev smoothness beta.
    """
    
    coeffs = [1]
    fns = [lambda x: 1]
    total = 0
    curr = 2
    Q = L**2/np.pi**(2*beta)
    while total <= Q:
        new_coeff = np.random.normal(0, 0.1)
        coeffs.append(new_coeff)
        if curr % 2 == 0:
            fns.append(lambda x: np.sqrt(2) * np.cos(np.pi*curr*x))
            total += new_coeff**2 * curr**(2*beta)
        else:
            fns.append(lambda x: np.sqrt(2) * np.sin(np.pi*(curr-1)*x))
            total += new_coeff**2 * (curr-1)**(2*beta)
        curr += 1
    return [coeffs, fns]

def linear_combination(coeffs, fns, x):
    return np.sum([coeffs[i]*fns[i](x) for i in range(len(coeffs))])
    
def uni_rejection_sample(coeffs, fns, n):
    """
    Rejection sampler for f(x) = \sum_i coeffs_i fns_i(x)
    Using uniform proposal distribution (which we hope will dominate 1/2 f(x))
    """
    x = []
    while len(x) < n:
        proposal_samples = np.random.uniform(0,1,n)
        us = np.random.uniform(0, 1, n)

        to_retain = [us[i] < linear_combination(coeffs, fns, proposal_samples[i])/2 for i in range(len(proposal_samples))]
        x.extend([proposal_samples[i] for i in range(len(proposal_samples)) if to_retain[i]])
    return x[0:n]

def rejection_sample(coeffs, fns, n, d):
    """
    Rejection sampler for f(x) = \sum_i coeffs_i fns_i(x)
    Using uniform proposal distribution (which we hope will dominate 1/2 f(x))
    """
    x = []
    while len(x) < n:
        proposal = np.random.uniform(0, 1, [n,d])
        us = np.random.uniform(0,1,n)
        to_retain = [i for i in range(n) if us[i] < linear_combination(coeffs, fns, [proposal[i,0], proposal[i,1]])/2]
        if x == []:
            x = proposal[to_retain,:]
        else:
            x = np.append(x, proposal[to_retain,:], 0)
    return x[0:n,:]

def plot_surface(coeffs, fns):
    """
    For a 2-d density, generate a mesh-grid of it's contours.
    """
    z = []
    x = np.arange(0, 1, 0.01)
    y = np.arange(0, 1, 0.01)
    for i in x:
        for j in y:
            z.append(linear_combination(coeffs, fns, [i,j]))
    z = np.array(z)
    Z = np.real(z.reshape(len(x), len(y)))
    X,Y = np.meshgrid(x,y)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X,Y,Z, cmap=cm.cool)
    plt.show()
    return (X,Y,Z)

def plot_fn_histogram(coeffs, fns):
    x = np.arange(0, 1, 0.01)
    y = np.arange(0, 1, 0.01)
    z = []
    for i in x:
        for j in y:
            z.append(linear_combination(coeffs, fns, [i,j]))
    z = np.array(z)
    Z = np.real(z.reshape(len(x), len(y)))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(Z, extent=[0, 1.0, 0, 1.0], interpolation='nearest')
    plt.show()

def plot_data_histogram(data):
    assert data.shape[1] == 2
    H, xedges, yedges = np.histogram2d(data[:,0], data[:,1], bins=50)
    H.shape, xedges.shape, yedges.shape = ((50, 50), (51,), (51,))
    
    extent = [yedges[0], yedges[-1], xedges[-1], xedges[0]]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.imshow(H, extent=extent, interpolation='nearest')
    plt.show()



if __name__=='__main__':
    [coeffs, fns] = trig_density(4, 5, 2)
    plot_surface(coeffs, fns)
    data = rejection_sample(coeffs, fns, 100000, 2)
    plot_fn_histogram(coeffs, fns)
    plot_data_histogram(data)
