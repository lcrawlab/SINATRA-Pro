#!/bin/python3

import numpy as np, sys
from scipy.special import expit
from scipy.integrate import trapz
from numpy.linalg import cholesky, solve
from scipy.stats import norm
from RATE import *

#from numba import jit
#@jit(nopython=True,parallel=True)
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

# Binary Laplace Apprxoimation classifier
# RW algorithm 3.1
def Laplace_Approximation(K,y,step_size=1.0):
    n = len(y)
    f = np.zeros(n,dtype=float)
    oldobj = 1e99
    for k in range(100):
        sigf = expit(f)
        W = np.diag(sigf*(1-sigf))
        sqrtW = np.sqrt(W)
        B = np.identity(n) + sqrtW @ K @ sqrtW
        L = cholesky(B)
        b = W @ f + y - sigf
        # Newton descent step
        a = b - solve(sqrtW @ L.T,solve(L,sqrtW @ K @ b))*step_size
        f = K @ a
        obj =  - .5 * (a.T @ f) - sigf
        newobj = np.mean(obj)
        print(newobj)        
        if np.fabs(newobj - oldobj) < 1e-5:
            break
        oldobj = newobj
    v = solve(L,sqrtW @ K)
    mean = f
    sigma = K - np.dot(v,v)
    log_marginal_likelihood = newobj - np.sum(np.log(np.diagonal(L)))
    return mean, sigma

def CovarianceFunction(K,sigma_n,X,y,xs,sigma=1.0):
    n = len(X)
    bandwidth = .5/sigma/sigma
    alpha = np.linalg.inv(K+np.diag(sigma_n)) @ y
    fx = 0
    kx = np.zeros(n)
    for i in range(n):
        kx[i] = np.exp(-np.mean((X[i]-xs)**2)*bandwidth)
        fx += alpha[i]*kx[i]
    return fx, kx

""" 
def Laplace_Approximation_Prediction(f,sigma_n,X,y,K,xs):
    n = len(y)
    sigf = expit(f)
    W = np.diag(sigf*(1-sigf))
    sqrtW = np.sqrt(W)
    B = np.identity(n) + sqrtW @ K @ sqrtW
    L = cholesky(B)
    fx, kx = CovarianceFunction(K,sigma_n,X,y,xs,1.0)
    fs = np.dot(kx,(y - sigf))
    v = solve(L, sqrtW @ kx)
    sigma = np.dot(kx,kx) - np.dot(v,v)
    z = np.linspace(fs-sigma*5,fs+sigma*5,201)
    predict_prob = trapz(expit(z)*norm.pdf(z,fs))/len(z)
    return predict_prob
"""

def Expectation_Propagation(K,y):
    n = len(y)
    mu = np.zeros(n,dtype=float)
    nu_tilde = np.zeros(n,dtype=float)
    tau_tilde = np.zeros(n,dtype=float)
    sigma = K
    for j in range(2):
        for i in range(n):
            a = sigma[i,i]**(-2.)

            tau_minus_i = a - tau_tilde[i]
            nu_minus_i = a*mu[i] - nu_tilde[i]
            mu_minus_i = nu_minus_i/tau_minus_i
            sigma_minus_i = tau_minus_i**(-.5)
            sigma_minus_i_sq = sigma_minus_i**2
             
            # Compute Marginal Moments
            z_i = y[i]*mu_minus_i/(np.sqrt(1+sigma_minus_i_sq))
            dnorm_z_i = norm.pdf(z_i)
            pnorm_z_i = norm.logcdf(z_i)

            mu_hat_i = mu_minus_i + (y[i]*sigma_minus_i_sq*dnorm_z_i)/(pnorm_z_i*np.sqrt(1+sigma_minus_i_sq))
            sigma_hat_i = np.sqrt(sigma_minus_i_sq-(sigma_minus_i_sq**2*dnorm_z_i)/((1+sigma_minus_i_sq)*pnorm_z_i)*(z_i + dnorm_z_i/(z_i)) )

            # Update Site Parameters
            delta_tau_tilde = sigma_hat_i**(-2) - tau_minus_i - tau_tilde[i]
            tau_tilde[i] = tau_tilde[i] + delta_tau_tilde
            nu_tilde[i] = sigma_hat_i**(-2)*mu_hat_i-nu_minus_i

            # Update Sigma, mu - the parameters of the posterior
            sigma = sigma - ((delta_tau_tilde**(-1)+sigma[i,i])**(-1))*(sigma[:,i].T @ sigma[:,i]) 
            mu = sigma @ nu_tilde

    #Recompute posterior parameters
    S_tilde = np.diag(tau_tilde)
    L = cholesky(np.identity(n) + np.sqrt(S_tilde) @ K @ np.sqrt(S_tilde))
    V = solve(L.T,np.sqrt(S_tilde) @ K)
    sigma = K - V.T @ V
    mu = sigma @ nu_tilde

    return mu, sigma

log_sqrt_2pi = .5*np.log(2*np.pi)
def log_likelihood(f,y):
    return np.sum(-(f*y)**2)-log_sqrt_2pi*len(f)

## Adopted from FastGP::ess
def Elliptical_Slice_Sampling(K,y,N_mcmc=100000,burn_in=1000,probit=True,seed=None):
    print("Running elliptical slice sampling...")
    if probit:
        log_lik = probit_log_likelihood
    else:
        log_lik = logistic_log_likelihood
    n = K.shape[0]
    N = y.size
    if isinstance(seed, int):
        np.random.seed(seed)
    mcmc_samples = np.zeros((burn_in+N_mcmc,N),dtype=float)
    norm_samples = np.random.multivariate_normal(mean = np.zeros(n), cov = K, size = burn_in+N_mcmc)
    unif_samples = np.random.uniform(low = 0, high = 1, size = burn_in+N_mcmc)
    theta = np.random.uniform(low = 0, high = 2*np.pi, size = burn_in+N_mcmc)
    theta_min = theta - 2*np.pi
    theta_max = theta + 2*np.pi
    for i in range(1,burn_in+N_mcmc):
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
    sys.stdout.write('\n')
    return mcmc_samples[burn_in:,:]

from time import perf_counter
def find_rate_variables_with_other_sampling_methods(X,y,bandwidth = 0.01,sampling_method = 'ESS',size = 100000,N_mcmc = 100000,burn_in = 1000,probit = True,seed = None, parallel = False, n_core = -1, verbose = False):
    n = X.shape[0]
    f = np.zeros(n)
    if verbose:
        sys.stdout.write('Calculating Covariance Matrix...\n')
    Kn = CovarianceMatrix(X.T,bandwidth)
    if sampling_method == 'Laplace':
        mu, sigma = Laplace_Approximation(Kn,y)
        samples = numpy.random.multivariate_normal(mean=mu,cov=sigma,size=size).T 
    elif sampling_method == 'EP':
        mu, sigma = Expectation_Propagation(Kn,y)
        samples = numpy.random.multivariate_normal(mean=mu,cov=sigma,size=size).T 
    elif sampling_method == 'ESS':
        samples = Elliptical_Slice_Sampling(Kn,y,N_mcmc=N_mcmc,burn_in=burn_in,probit=probit,seed=seed)
    else:
        return
    kld, rates, delta, eff_samp_size = RATE(X=X,f_draws=samples,parallel=parallel,n_core=n_core,verbose=verbose)
    return kld, rates, delta, eff_samp_size
 
