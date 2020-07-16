#!/bin/python3

#from mesh import *
#from euler import *
#from directions import *
from gp import *

#import matplotlib.pyplot as plt

import numpy as np
import os


from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF

X = []
y = []
for i in range(1,91):
    ec = np.loadtxt('ec_sphere/WT/WT_chimera_%d.dat'%i)[1:].flatten()
    if ec[-1] > 100 and ec[-1] < 200:
        X.append(ec)
        y.append(0)
for i in range(1,91):
    ec = np.loadtxt('ec_sphere/R164S/R164S_chimera_%d.dat'%i)[1:].flatten()
    if ec[-1] > 100 and ec[-1] < 200:
        X.append(ec)
        y.append(1)
#for i in range(1,91):
#    ec = np.loadtxt('ec_sphere_2/R164G/R164G_chimera_%d.dat'%i)[1:].flatten()
#    if ec[-1] > 100 and ec[-1] < 200:
#        X.append(ec)
#        y.append(2)

X = np.array(X)
y = np.array(y)

Kn = GuassianKernel(X,bandwidth=0.02)
mean,covariance = LaplaceApproximation(Kn,y)
#print(mean,len(mean))
#print(covariance,covariance.shape)

'''
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





