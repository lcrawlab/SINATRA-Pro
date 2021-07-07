#!/bin/python3

import numpy as np, sys
from scipy.stats import norm
from RATE import *

def CovarianceMatrix(x,bandwidth=0.01):
    bandwidth = 1./(2*bandwidth**2)
    n = x.shape[1]
    K = np.zeros((n,n),dtype=float)
    for i in range(n):
        K[i,i] = 1
        for j in range(i+1,n):
            K[i,j] = np.exp(-np.mean((x[:,i]-x[:,j])**2)*bandwidth)
            K[j,i] = K[i,j]
    return K

def probit_log_likelihood(latent_variables, class_labels):
    return np.sum(np.log(norm.cdf(latent_variables*class_labels)))

def logistic_log_likelihood(latent_variables, class_labels):
    return(-np.sum(np.log(1.0+np.exp(latent_variables*class_labels))))

## Adopted from FastGP::ess
def Elliptical_Slice_Sampling(K,y,n_mcmc=100000,burn_in=1000,probit=True,seed=None,verbose=False):
    if verbose:
        print("Running elliptical slice sampling...")
    if probit:
        log_lik = probit_log_likelihood
    else:
        log_lik = logistic_log_likelihood
    n = K.shape[0]
    N = y.size
    if isinstance(seed, int):
        np.random.seed(seed)
    mcmc_samples = np.zeros((burn_in+n_mcmc,N),dtype=float)
    norm_samples = np.random.multivariate_normal(mean = np.zeros(n), cov = K, size = burn_in+n_mcmc)
    unif_samples = np.random.uniform(low = 0, high = 1, size = burn_in+n_mcmc)
    theta = np.random.uniform(low = 0, high = 2*np.pi, size = burn_in+n_mcmc)
    theta_min = theta - 2*np.pi
    theta_max = theta + 2*np.pi
    for i in range(1,burn_in+n_mcmc):
        if verbose:
            if i < burn_in:
                sys.stdout.write('Burning in...\r')
            else:
                sys.stdout.write('Elliptical slice sampling Step %d...\r'%(i-burn_in+1))
            sys.stdout.flush()
        f = mcmc_samples[i-1,:]
        llh_thresh = log_lik(f,y) + np.log(unif_samples[i])
        f_star = f * np.cos(theta[i]) + norm_samples[i,:] * np.sin(theta[i])
        while(log_lik(f_star,y) < llh_thresh):
            if theta[i] < 0:
                theta_min[i] = theta[i]
            else:
                theta_max[i] = theta[i]
            theta[i] = np.random.uniform(low = theta_min[i],high = theta_max[i], size = 1)
            f_star = f * np.cos(theta[i]) + norm_samples[i,:] * np.sin(theta[i])
        mcmc_samples[i,:] = f_star
    if verbose:
        sys.stdout.write('\n')
    return mcmc_samples[burn_in:,:]

def find_rate_variables_with_other_sampling_methods(X,y,bandwidth = 0.01,sampling_method = 'ESS', n_mcmc = 100000,burn_in = 1000,probit = True,seed = None, parallel = False, n_core = -1, verbose = False):
    n = X.shape[0]
    f = np.zeros(n)
    if verbose:
        sys.stdout.write('Calculating Covariance Matrix...\n')
    Kn = CovarianceMatrix(X.T,bandwidth)
    samples = Elliptical_Slice_Sampling(Kn,y,n_mcmc=n_mcmc,burn_in=burn_in,probit=probit,seed=seed,verbose=verbose)
    kld, rates, delta, eff_samp_size = RATE(X=X,f_draws=samples,parallel=parallel,n_core=n_core,verbose=verbose)
    return kld, rates, delta, eff_samp_size
 
