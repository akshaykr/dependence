from scipy.integrate import dblquad, quad

def numeric_integration(f, l_limit, u_limit):
    """
    A wrapper function for python's quad and dblquad. f is the function handle.
    l_limit and u_limit are lists/ tuples which specify the lower and upper
    limits along each limit.
    """
    num_dims = len(l_limit);
    if num_dims == 1:
      ret = quad(f, l_limit[0], u_limit[0]);
    elif num_dims == 2:
      ret = dblquad(f, l_limit[0], u_limit[0], lambda x: l_limit[1],
                    lambda x: u_limit[1]);
    return ret[0];

