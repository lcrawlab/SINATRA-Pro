#!/bin/python3

import numpy as np
from scipy.special import expit
from numpy.linalg import cholesky, solve
from scipy.stats import norm
#import time

from numba import jit

@jit(nopython=True,parallel=True)
def GuassianKernel(x,bandwidth=1.0):
    n = x.shape[0]
    K = np.zeros((n,n))
    for i in range(n):
        K[i,i] = 1
        for j in range(i+1,n):
            y = np.mean((x[i]-x[j])**2)*bandwidth
            K[i,j] = np.exp(-y)
            K[j,i] = K[i,j]
    return K

# Binary Laplace Apprxoimation classifier
# RW algorithm 3.1
def LaplaceApproximation(K,y,step_size=0.1):
    n = len(y)
    f = np.zeros(n,dtype=float)
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
        obj = - .5 * (a.T @ f) - sigf
        print(np.mean(obj))
    v = solve(L,sqrtW @ K)
    mean = f
    covariance = K - np.dot(v,v)
    return mean, covariance


def ExpectationPropagation(K,y):
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
            print(z_i)
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

