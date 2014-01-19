import numpy as np
from numpy import power as powr

def leg_poly_0(u):
  return 1/np.sqrt(2) * np.ones(np.shape(u));

def leg_poly_1(u):
  return np.sqrt(3/2) * u;

def leg_poly_2(u):
  return 1/2.0 * np.sqrt(5/2) * ( 3*powr(u,2) - 1);

def leg_poly_3(u):
  return 1/2.0 * np.sqrt(7/2) * (5*powr(u,3) - 3*u);

def leg_poly_4(u):
  return 1/(128*np.sqrt(2)) * (384*powr(u,4) + 
                               1152*np.multiply((powr(u,2) - 1), powr(u,2)) + 
                               144*powr((powr(u, 2) - 1), 2) );

def leg_poly_5(u):
  return 1/3840.0 * np.sqrt(11/2) * (3840*powr(u,5) +
                      7200 * np.multiply(powr((powr(u,2) - 1), 2), u) + 
                      19200 * np.multiply((powr(u,2) - 1), powr(u,3)) );

def kernel_1D(x, c, h, order):
  basis_fn = { 0: leg_poly_0,
               1: leg_poly_1,
               2: leg_poly_2,
               3: leg_poly_3,
               4: leg_poly_4,
               5: leg_poly_5
             };
  u = (x - c)/h;
  ret = 0.0;
  for i in range(0, order+1):
    ret = ret + basis_fn[i](0) * basis_fn[i](u);
  ret = np.multiply(ret, (abs(u) < 1).astype(float));
  return ret;

def kernel(X, order, h=1.0, centre=None):
  """
  X is an num_data x num_dims data arra.
  order is the order of the kernel. (1 <= order <= 5)
  h is the bandwidth
  centre is a 1 x num_dims array.
  The function returns a num_data x 1 array where the ith row is
  K((x[i,:] - c)/h)
  """
  # Prelims
  order = int(order);
  h = float(h);
  X = np.matrix(X);
  (num_data, num_dims) = np.shape(np.matrix(X));
  if centre is None:
    centre = np.zeros((1, num_dims));

  ret = np.ones((num_data, 1));
  for dim_idx in range(0, num_dims):
    ret *= kernel_1D(X[:, dim_idx], centre[0, dim_idx], h, order);
  return ret 

