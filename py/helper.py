from scipy.integrate import dblquad, quad
import numpy as np

"""
Helper routines for numeric integration.
"""

def numeric_integration(f, l_limit, u_limit):
    """
    A wrapper function for python's quad and dblquad. f is the function handle.
    l_limit and u_limit are lists/ tuples which specify the lower and upper
    limits along each limit.
    """
    num_dims = len(l_limit);
    if num_dims == 1:
      ret = quad(f, l_limit[0], u_limit[0], full_output=1);
    elif num_dims == 2:
      ret = dblquad(f, l_limit[0], u_limit[0], lambda x: l_limit[1],
                    lambda x: u_limit[1]);
    return ret[0];


def fast_integration(f, l_limit, u_limit, pts = 1000):
    """
    Faster integration by splitting the domain into rectangles.

    f is the function handle. In contrast with the numeric_integration routine, this is a function of one variable which is a n-by-d matrix of points to evaluate.
    l_limit is a vector of the lower limits
    u_limit is the vector of the upper limits
    """
    num_dims = len(l_limit)
    if num_dims == 1:
        pts = np.matrix(np.arange(l_limit[0], u_limit[0], (u_limit[0]-l_limit[0])/pts)).T
        return np.mean(f(pts))

    if num_dims == 2:
        pts0 = np.arange(l_limit[0], u_limit[0], 0.05)
        pts1 = np.arange(l_limit[1], u_limit[1], 0.05)
        X,Y = np.meshgrid(pts0, pts1)

        v = np.matrix(zip(X.reshape(len(pts0)*len(pts1),), Y.reshape(len(pts0)*len(pts1),)))

        return np.mean(f(v))
