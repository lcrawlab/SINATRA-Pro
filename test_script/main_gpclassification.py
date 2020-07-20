#!/bin/python3

#from mesh import *
#from euler import *
#from directions import *
from gp import *

#import matplotlib.pyplot as plt

import numpy as np
import os

X = []
y = []

prefixes = ['ec_sphere/WT/WT_chimera_','ec_sphere/R164S/R164S_chimera_']#,'ec_sphere/R164G/R164G_chimera_']
n_sample = 50
class_labels = [1,-1]

for i in range(len(prefixes)):
    for j in range(1,n_sample+1):
        ec = np.loadtxt('%s%d.dat'%(prefixes[i],j))[1:].flatten()
        if ec[-1] > 100 and ec[-1] < 200:
            X.append(ec)
            y.append(class_labels[i])
X = np.array(X)
y = np.array(y)

# Hard coded Laplace Approximation
K = GuassianKernel(X,sigma=1.0)
mean, sigma_n, f, log_marginal_likelihood = LaplaceApproximation(K,y,step_size=0.001)
#print(mean, sigma, f, log_marginal_likelihood)

xs = np.loadtxt('ec_sphere/WT/WT_chimera_100.dat')[1:].flatten()
LaplaceApproximationPrediction(f,sigma_n,X,y,K,xs)

#mean, sigma = ExpectationPropagation(K,y)
#print(mean,len(mean))
#print(sigma,sigma.shape)


'''
import pymc3 as pm
latent_gp_model = pm.Model()

conv_func = pm.gp.cov.ExpQuad(1, ls=0.1)
gp = pm.gp.Latent(conv_func=conv_func)
f = gp.prior("f",X=X)
'''


'''
## sk learn GPC

from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF

kernel = 1.0*RBF(1.0)
gpc = GaussianProcessClassifier(kernel=kernel,random_state=0,max_iter_predict=1000).fit(X,y)

gpcscore = gpc.score(X, y)
print(gpcscore)
X = []
for i in range(1,101):
    ec = np.loadtxt('ec_sphere_2/WT/WT_chimera_%d.dat'%i)[1:].flatten()
    X.append(ec)
X = np.array(X)

print(gpc.predict_proba(X[:,:]))
print(gpc.predict(X[:,:]))
'''


