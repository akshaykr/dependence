import numpy as np

"""
Module for generating densities and related subroutines.
"""

def trigonometric_uni_density(beta, L = 10):
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
#         if curr == 1:
#             fns.append(lambda x: 1)
#             total += new_coeff**2 * curr**(2*beta)
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
    
def rejection_sample(coeffs, fns, n):
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
