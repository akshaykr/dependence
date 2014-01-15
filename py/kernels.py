import numpy as np


def kernel_1_1D(x, h=1.0, c=0.0):
  """
  Implements a 1D first order Legendre Kernel.
  x is the point at which the kernel should be evaluated. h is the bandwith.
  c is the centre. The function returns (1/h) K((x-c)/h).
  """
  u = (x - c)/h;
  return (1/h) * (1/2.0) * float(abs(u) < 1);

def kernel_2_1D(x, h=1.0, c =0.0):
  """
  Implements a 1D second order Legendre Kernel.
  Inputs/ Outputs same as kernel_1_1D.
  """
  u = (x - c)/h;
  return (1/h) * (1/8.0) * (9 - 15*u**2) * float(abs(u) < 1);
  
def kernel_1_2D(x, h=1.0, c=[0.0, 0.0]):
  """
  Implements a 2D first order Legendre Kernel.
  x is the point at which the kernel should be evaluated. h is the bandwith.
  c is the centre. The function returns (1/h**2) K(||x-c||/h).
  """
  return kernel_1_1D(x[0], h, c[0]) * kernel_1_1D(x[1], h, c[1]);

def kernel_2_2D(x, h=1.0, c=[0.0, 0.0]):
  """
  Implements a 2D second order Legendre Kernel.
  Inputs/ Outputs same as kernel_1_2D
  """
  return kernel_2_1D(x[0], h, c[0]) * kernel_2_1D(x[1], h, c[1]);
  
def kernel(x, order, h=1.0, c=None):
  """
  Just a wrapper function. Chooses one of the above 4 kernels depending on dim
  and order.
  """
  if c is None:
    c = np.zeros(len(x));
  options = { (1, 1): kernel_1_1D, (1,2): kernel_1_2D, 
              (2, 1): kernel_2_1D, (2,2): kernel_2_2D };
  return options[(order, len(x))](x, h, c);

