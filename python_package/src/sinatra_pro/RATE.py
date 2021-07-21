#!/bin/python3

import sys
import numpy as np
from scipy.linalg import pinv

def sherman_r(A,u,v):
    """ Sherman-Morrisoni formula to compute the inverse of the sum of an invertible matrix A and the outer product of vectors u and v."""
    x = v.T @ A @ u + 1
    return A - ((A @ u) @ (v.T @ A)) * (1./x)

def calc_kld(mu,Lambda,V,q,verbose=False):
    """Calculate KLD for column q """
    if verbose:
        sys.stdout.write("Calculating KLD(%d)...\r"%q)
        sys.stdout.flush()
    U_Lambda_sub = sherman_r(Lambda,V[:,q],V[:,q].T)
    U_no_q = np.delete(U_Lambda_sub,q,0)
    U_no_qq = np.delete(U_no_q,q,1)
    alpha = U_no_q[:,q].T @ U_no_qq @ U_no_q[:,q]
    kld = mu[q]**2 * alpha * .5
    return kld

def RATE(X,f_draws=None,prop_var=1,low_rank=False,parallel=False,n_core=-1,verbose=False):    
    """
    Variable Prioritization via RelATive cEntrality (RATE) centrality measures.

    This function assumes that one has already obtained (posterior) draws/estimates of a nonparametric or nonlinear function as suggested in Crawford et al. (2018)

    'X' is the nxp design matrix (e.g. genotypes) where n is the number of samples and p is the number of dimensions. This is the original input data
    
    'f_draws' is the Bxn matrix of the nonparametric model estimates (i.e. f.hat) with B being the number of sampled (posterior) draws;
    
    'prop_var' is the desired proportion of variance that the user wants to explain when applying singular value decomposition (SVD) to the design matrix X (this is preset to 1);
    
    'low_rank' is a boolean variable detailing if the function will use low rank matrix approximations to compute the RATE values --- note that this highly recommended in the case that the number of covariates (e.g. SNPs, genetic markers) is large; 
       
    If `parallel` is set to True, the program runs on multiple cores for RATE calculations, 
    then `n_core` will be the number of cores used (the program uses all detected cores if `n_core` is not provided).

    If `verbose` is set to True, the program prints progress on command prompt.
 
    """
    if verbose:
        sys.stdout.write("Calculating RATE...\n")

    if parallel:
        import multiprocessing
        from joblib import Parallel, delayed
        if n_core == -1:
            n_core = multiprocessing.cpu_count()

    if low_rank:
        ### Take the SVD of the Design Matrix for Low Rank Approximation ###
        u, s, vh = np.linalg.svd(X,full_matrices=False,compute_uv=True)
        dx = s > 1e-10
        s_sq = s**2
        px = np.cumsum(s_sq/np.sum(s_sq)) < prop_var
        r_X = np.logical_and(dx,px)
        u = ((1. / s[r_X]) * u[:,r_X]).T
        v = vh.T[:,r_X]
        # Now, calculate Sigma_star
        SigmaFhat = np.cov(f_draws, rowvar=False)
        Sigma_star = u @ SigmaFhat @ u.T 
        # Now, calculate U st Lambda = U %*% t(U)
        u_Sigma_star, s_Sigma_star, vh_Sigma_star = np.linalg.svd(Sigma_star,full_matrices=False,compute_uv=True)
        r = s_Sigma_star > 1e-10
        tmp = 1./np.sqrt(s_Sigma_star[r]) * u_Sigma_star[:,r].T
        U = pinv(v).T @ tmp.T
        V = v @ Sigma_star @ v.T #Variances
        mu = v @ u @ np.average(f_draws,axis=0) #Effect Size Analogues
    else:
        beta_draws = (pinv(X) @ f_draws.T).T
        V = np.cov(beta_draws,rowvar=False)
        D = pinv(V)
        D_u, D_s, D_vh = np.linalg.svd(D,full_matrices=False,compute_uv=True)
        r = np.sum(D_s > 1e-10)
        U = np.multiply(np.sqrt(D_s[:r]),D_u[:,:r])
        mu = np.average(beta_draws,axis=0)
    
    mu = np.fabs(mu)

    ### Create Lambda ###
    Lambda = U @ U.T

    ### Compute the Kullback-Leibler divergence (KLD) for Each Predictor ###
    if parallel:
        kld = Parallel(n_jobs=n_core)(delayed(calc_kld)(mu,Lambda,V,q,verbose) for q in range(mu.size))     
        kld = np.array(kld)
    else:
        kld = np.zeros(mu.size,dtype=float)
        for q in range(mu.size):
            kld[q] = calc_kld(mu,Lambda,V,q,verbose)
    if verbose: 
        sys.stdout.write("\n")
        sys.stdout.write("KLD calculation Completed.\n")
    
    ### Compute the corresponding “RelATive cEntrality” (RATE) measure ###
    rates = kld/np.sum(kld)

    ### Find the entropic deviation from a uniform distribution ###
    #delta = np.sum(rates*np.log((len(mu)-len(nullify))*rates))
    delta = np.sum(rates*np.log(len(mu)*rates))

    ### Calibrate Delta via the effective sample size (ESS) measures from importance sampling ###
    #(Gruber and West, 2016, 2017)
    eff_samp_size = 1./(1.+delta)*100.

    ### Return a list of the values and results ###
    return kld, rates, delta, eff_samp_size


