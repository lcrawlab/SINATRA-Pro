#!/bin/python3

import numpy as np
from scipy.special import expit
from numpy.linalg import solve

#import time

from numba import jit

@jit(nopython=True,parallel=True)
def GuassianKernel(x,bandwidth=1.0):
    n = x.shape[0]
    K = np.zeros((n,n),dtype=float)
    for i in range(n):
        K[i,i] = 1
        for j in range(i+1,n):
            y = np.mean((x[i]-x[j])**2)*bandwidth
            K[i,j] = np.exp(-y)
            K[j,i] = K[i,j]
    return K

def LaplaceApproximation(Kn,y):
    n = len(y)
    f = np.zeros(n,dtype=float)
    oldobj = 0
    for k in range(100):
        #print(k)
        sigf = expit(f)
        W = np.diag(sigf*(1-sigf))
        sqrtW = np.sqrt(W)
        B = np.identity(n) + sqrtW @ Kn @ sqrtW
        L = np.linalg.cholesky(B)
        b = W @ f + y - sigf
        a = b - solve(sqrtW @ L.T,solve(L,sqrtW @ Kn @ b))*0.1
        f = Kn @ a
        obj = - .5 * (a.T @ f) - sigf
        print(np.mean(obj-oldobj))
        oldobj = obj
    v = solve(L,sqrtW @ Kn)
    mean = f
    covariance = Kn - np.dot(v,v)
    return mean, covariance



