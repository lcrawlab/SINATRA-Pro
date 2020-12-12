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

def RATE(X,f_draws=None,pre_specify=False,beta_draws=None,prop_var=1,nullify=[],snp_nms=None,n_core=1):
    if n_core == -1:
        n_core = multiprocessing.cpu_count()
    ### Take the SVD of the Design Matrix for Low Rank Approximation ###
    u, s, vh = np.linalg.svd(X,full_matrices=False,compute_uv=True)
    dx = s > 1e-10
    s_sq = s**2
    px = np.cumsum(s_sq/np.sum(s_sq)) < prop_var
    r_X = np.logical_and(dx,px)
    u = (1. / s[r_X] * u[:,r_X]).T
    v = vh.T[:,r_X]
    
    # Now, calculate Sigma_star
    SigmaFhat = np.cov(f_draws, rowvar=False)
    Sigma_star = u @ SigmaFhat @ u.T
     
    # Now, calculate U st Lambda = U %*% t(U)
    u_Sigma_star, s_Sigma_star, vh_Sigma_star = np.linalg.svd(Sigma_star,full_matrices=False,compute_uv=True)
    r = s_Sigma_star > 1e-10

    tmp = 1./np.sqrt(s_Sigma_star[r]) * u_Sigma_star[:,r].T
    U = np.linalg.pinv(v).T @ tmp.T
    
    V = v @ Sigma_star @ v.T #Variances
    mu = v @ u @ np.average(f_draws,axis=0) #Effect Size Analogues
    ### Create Lambda ###
    Lambda = U @ U.T
        
    ### Compute the Kullback-Leibler divergence (KLD) for Each Predictor ###
    kld = np.zeros(mu.size,dtype=float)
    for i in range(mu.size):
        if i in nullify:
            continue
        q = np.unique(np.append(i,nullify)).astype(int)
        m = np.fabs(mu[q])
        U_Lambda_sub = woodbury(Lambda,V[:,q],V[:,q].T)
        U_no_q = np.delete(U_Lambda_sub,q,0)
        U_no_qq = np.delete(U_no_q,q,1)
        alpha = U_no_q[:,q].T @ U_no_qq @ U_no_q[:,q]
        kld[i] = (m.T @ alpha @ m) * .5
    
    ### Compute the corresponding “RelATive cEntrality” (RATE) measure ###
    rates = kld/np.sum(kld)

    ### Find the entropic deviation from a uniform distribution ###
    delta = np.sum(rates*np.log((len(mu)-len(nullify))*rates))

    ### Calibrate Delta via the effective sample size (ESS) measures from importance sampling ###
    #(Gruber and West, 2016, 2017)
    eff_samp_size = 1./(1.+delta)*100.

    ### Return a list of the values and results ###
    return kld, rates, delta, eff_samp_size


