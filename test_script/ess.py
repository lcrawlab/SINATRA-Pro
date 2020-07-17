#!/bin/python3

import numpy as np
from scipy.stats import norm

log_sqrt_2pi = np.log(2*np.pi)/2.
def log_likelihood(f, y):
    return np.sum(-(f*y)**2/2.+log_sqrt_2pi)

def EllipticalSliceSampling(K,y):
    

