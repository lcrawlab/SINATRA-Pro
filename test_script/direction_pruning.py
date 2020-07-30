#/bin/python3

import numpy as np

def DataCorrelation(X1,X2):
    n = len(X1)
    m = len(X2)
    cor = np.zeros((n,m))
    for i in range(n):
        for j in range(m):
            cor[i,j] = np.correlate(X1[i],X2[j])
            #cor[i,j] = np.corrcoef(X1[i],X2[j])
    return cor

def DirectionCorrelation(X,n_direction):
    n_sample = X.shape[0]
    n_filtration = int(X.shape[1]/n_direction)
    corr = np.zeros((n_direction,n_direction))
    for i in range(n_direction):
        for j in range(i,n_direction):
            x1 = np.mean(X[:,i*n_filtration:(i+1)*n_filtration],axis=0)
            x2 = np.mean(X[:,j*n_filtration:(j+1)*n_filtration],axis=0)
            corr[i,j] = np.correlate(x1,x2)
    return corr