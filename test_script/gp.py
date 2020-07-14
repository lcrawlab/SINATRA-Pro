#!/bin/python3

import numpy as np
from scipy.spatial.distance import pdist, squareform

import time

def GuassianKernal(x,bandwidth):
    n = x.shape[0]
    K = np.zeros((n,n),dtype=float)
    for i in range(n):
        K[i,i] = 1
        for j in range(i+1,n):
            y = np.mean((x[i]-x[j])**2)*bandwidth
            K[i,j] = np.exp(-y)
            K[j,i] = K[i,j]
    return K


x = np.arange(2000).reshape(100,20)
start = time.time()
K = GuassianKernal(x,0.1)
end = time.time()
print(end-start)

print(K)





