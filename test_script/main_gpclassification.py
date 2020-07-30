#!/bin/python3

#from mesh import *
#from euler import *
#from directions import *
from gp import *
from direction_pruning import *

#import matplotlib.pyplot as plt

import numpy as np
import os


mutants = ['WT','R164G','R164S']
#mutants = ['WT']

'''
for i in range(len(mutants)):
    inliers = []
    for j in range(1,101):
        ec = np.loadtxt('ec_sphere/%s/%s_chimera_%d.dat'%(mutants[i],mutants[i],j))[1:].flatten()
        if ec [-1] > 50 and ec[-1] < 200:
            inliers.append(j-1)
    np.savetxt('inliers_%s.txt'%mutants[i],inliers)
'''

X = []
y = []
class_labels = [1,-1]
for i in range(len(mutants)):
    #inliers = np.loadtxt('inliers_%s.txt'%mutants[i])
    for j in range(100):
        ec = np.loadtxt('dect_sphere/%s/%s_chimera_%d.dat'%(mutants[i],mutants[i],j+1))[1:].flatten()
        X.append(ec)
        #y.append(class_labels[i])

X = np.array(X)
#y = np.array(y)


corr = np.corrcoef(X)
#corr = [a[np.isfinite(a)] for a in corr]
#corr = [a for a in corr if len(a) > 0]

from matplotlib import pyplot as plt
from matplotlib.cm import get_cmap

fig = plt.figure(figsize=(8,6))

ax = fig.add_axes([0.1,0.1,0.8,0.8])
im = ax.imshow(corr,aspect=1.0,interpolation ='none',vmax=1,vmin=-1,cmap=get_cmap('seismic'))

cax = fig.add_axes([0.825, 0.1, 0.04, 0.8])
fig.colorbar(im, cax=cax, orientation='vertical')
cax.set_ylabel('Pearson\'s correlation coefficient',fontsize=12)
cax.tick_params(labelsize=10)
plt.savefig('correlation_mesh.pdf')
plt.show()
plt.close()


'''
# Hard coded Laplace Approximation
K = GuassianKernel(X,sigma=1.0)
mean, sigma_n, f, log_marginal_likelihood = LaplaceApproximation(K,y,step_size=0.01)
print(mean, sigma_n, f, log_marginal_likelihood)

xs = np.loadtxt('ec_sphere/R164S/R164S_chimera_98.dat')[1:].flatten()
predict_prob = LaplaceApproximationPrediction(f,sigma_n,X,y,K,xs)
print(predict_prob)

#mean, sigma = ExpectationPropagation(K,y)
#print(mean,len(mean))
#print(sigma,sigma.shape)
'''

'''
import pymc3 as pm

latent_gp_model = pm.Model()

conv_func = pm.gp.cov.ExpQuad(1, ls=0.1)
gp = pm.gp.Latent(mean_func, cov_func)
f = gp.prior("f",X=X)
f_star = gp.conditional("f_star", X_star)
'''
'''
# sklearn GPC

from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF

kernel = 1.0*RBF(1.0)
gpc = GaussianProcessClassifier(kernel=kernel,random_state=0,max_iter_predict=2000).fit(X,y)

gpcscore = gpc.score(X, y)
print(gpcscore)

print(gpc.predict_proba(X[:,:]))
print(gpc.predict(X[:,:]))
'''

