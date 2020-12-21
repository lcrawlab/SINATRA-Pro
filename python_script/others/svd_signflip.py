#!/bin/python3

import matplotlib.pyplot as plt
import numpy as np
import os, sys
from scipy.linalg import svd, pinv, pinv2

x = np.loadtxt("data/dect_WT_R164S_65_213_with_h_4_15_01_4_25_norm.txt")
pinvX = np.linalg.pinv(x)
print(pinvX)
print(pinv(x))
print(pinv2(x))
exit()

def SignFlip(X,u,s,vh):
    N = s.shape[0]
    sk_left = np.zeros(N,dtype=float)
    sk_right = np.zeros(N,dtype=float)
    for k in range(N):
        s_nok = np.copy(s)
        s_nok[k] = 0
        y = X - (u * s_nok) @ vh
        for j in range(y.shape[1]):
            uk_yj = u[k] @ y[:,j]
            sk_left[k] += np.sign(uk_yj)*(uk_yj)**2
        for i in range(y.shape[0]):
            vk_yi = vh[k] @ y[i,:]
            sk_right[k] += np.sign(vk_yi)*(vk_yi)**2
    for k in range(N):
        if sk_left[k] * sk_right[k] < 0:
            if sk_left[k] < sk_right[k]:
                sk_left[k] *= -1
            else:
                sk_right[k] *= -1
        u[k] = np.sign(sk_left[k])*u[k]
        vh[k] = np.sign(sk_right[k])*vh[k]
    return u, vh

u, s, vh = svd(x,full_matrices=False,compute_uv=True,lapack_driver='gesvd')
np.save("u.npy",u)
np.save("vh.npy",vh)
np.save("s.npy",s)

#u = np.load("u.npy")
#vh = np.load("vh.npy")
#s = np.load("s.npy")
#u_py, vh_py = SignFlip(x,u,s,vh)
#print(u, vh)

u_r = np.loadtxt("data/svd_u_R.dat")
vh_r = np.loadtxt("data/svd_v_R.dat").T
s_r = np.loadtxt("data/svd_d_R.dat")

#u, vh = SignFlip(x,u,s,vh)
print(u_r/u, vh_r/vh)

