#!/bin/python3

#from mesh import *
#from euler import *
#from directions import *
from gp import *

import matplotlib.pyplot as plt
import numpy as np
import os, sys

x = np.loadtxt("data/dect_WT_R164S_65_213_4_15_01_4_25_norm.txt")
N = x.shape[0]
y = np.zeros(N,dtype=int)
y[:int(N/2)].fill(-1)
y[int(N/2):].fill(1)

#rates = find_rate_variables_with_other_sampling_methods(x,y,bandwidth=0.01,sampling_method="ESS")

#Kn = CovarianceMatrix(x.T,bandwidth=0.01)
#np.savetxt('Kn.dat',Kn)
#Kn = np.loadtxt('Kn.dat')
#samples = Elliptical_Slice_Sampling(Kn,y,N_mcmc=100000,burn_in=1000,probit=True)
#np.savetxt('ess_samples.dat',samples)
#samples = np.loadtxt('ess_samples.dat')
#np.save('ess_samples.npy',samples)

samples = np.load('data/ess_samples.npy')
#print(samples.shape)
#mean = np.mean(samples,axis=0)
#std = np.std(samples,axis=0)
#print(mean,std)
#plt.plot(range(len(mean)),mean)
#plt.fill_between(range(len(mean)),mean-std,mean+std)
#plt.show()
kld, rates, delta, eff_samp_size = RATE(x,f_draws=samples)
np.savetxt('data/rates.txt',rates)
print(kld,rates,delta,eff_samp_size)


