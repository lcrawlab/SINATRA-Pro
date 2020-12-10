#!/bin/python3

import numpy as np
import time
import multiprocessing
from joblib import Parallel, delayed

# Woodbury matrix idenity
def woodbury(A, u, v):
    A_inv_diag = 1./np.diag(A)
    B_inv = np.linalg.inv((v * A_inv_diag) @ u + np.eye(u.shape[1]))
    return np.diag(A_inv_diag) - (A_inv_diag.reshape(-1,1) * u @ B_inv @ v * A_inv_diag)

"""
p = 10000
k = 100
A = np.diag(np.random.randn(p))
u = np.random.randn(p,k)
v = np.random.randn(k,p)

vbegin = time.perf_counter()
print(woodbury(A,u,v))
end = time.perf_counter()
print(end-begin)
"""

def RATE(X,f_draws=None,pre_specify=False,beta_draws=None,prop_var=1,nullify=None,snp_nms=None,cores=1):
    ### Take the SVD of the Design Matrix for Low Rank Approximation ###
    u, s, vh = np.linalg.svd(X,full_matrices=False,compute_uv=True)
    dx = s > 1e-10
    s_sq = s**2
    px = np.cumsum(s_sq/np.sum(s_sq)) < prop_var
    r_X = np.logical_and(dx,px)
    u = 1./s[r_X]*u[:,r_X]
    v = vh.T[:,r_X]
    
    # Now, calculate Sigma_star
    SigmaFhat = np.cov(f_draws)
    Sigma_star = u @ SigmaFhat @ u.T

    print(Sigma_star)
     
    # Now, calculate U st Lambda = U %*% t(U)
    u_Sigma_star, s_Sigma_star, vh_Sigma_star = np.linalg.svd(Sigma_star,full_matrices=False,compute_uv=True)
    r = s_Sigma_star > 1e-10
    
    tmp = 1./np.sqrt(s_Sigma_star[r])*u_Sigma_star[,r].T
    U = np.linalg.pinv(v).T @ tmp.T
    
    """
    V = v%*%Sigma_star%*%t(v) #Variances
    mu = v%*%u%*%colMeans(f.draws) #Effect Size Analogues
    ### Create Lambda ###
    Lambda = tcrossprod(U)
      
    ### Compute the Kullback-Leibler divergence (KLD) for Each Predictor ###
    int = 1:length(mu); l = nullify;
        
    if(length(l)>0){int = int[-l]}
        
    KLD = foreach(j = int, .combine='c')%dopar%{
         q = unique(c(j,l))
         m = abs(mu[q])
                                  
    U_Lambda_sub = sherman_r(Lambda,V[,q],V[,q])
    alpha = t(U_Lambda_sub[-q,q])%*%U_Lambda_sub[-q,-q]%*%U_Lambda_sub[-q,q]
    kld = (t(m)%*%alpha%*%m)/2
    names(kld) = snp.nms[j]
    return kld
    """

#X = np.loadtxt('../../ubq/WT/restrained_all/tda/dect_D_A_0_4_15_01_8_50_norm.txt')

