from scipy.integrate import dblquad, quad
import numpy as np

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


def fast_integration(f, l_limit, u_limit):
    num_dims = len(l_limit)
    if num_dims == 1:
        pts = np.matrix(np.arange(l_limit[0], u_limit[0], 0.01)).T
        return np.mean(f(pts))

    if num_dims == 2:
        pts0 = np.arange(l_limit[0], u_limit[0], 0.01)
        pts1 = np.arange(l_limit[1], u_limit[1], 0.01)
        X,Y = np.meshgrid(pts0, pts1)

        v = np.matrix(zip(X.reshape(len(pts0)*len(pts1),), Y.reshape(len(pts0)*len(pts1),)))

        return np.mean(f(v[:,0], v[:,1]))
