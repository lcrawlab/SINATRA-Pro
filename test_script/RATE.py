#!/bin/python3

import numpy as np

def RATE(X,f.draws=Null,pre.specify=False,beta.draws=Null,prop_var=1,nullify=NULL,snp.nms=NULL,cores=1)
    ### Take the SVD of the Design Matrix for Low Rank Approximation ###
    u, s, vh = numpy.linalg.svd(X)
    dx = s > 1e-10
    s_sq = s**2
    px = np.cumsum(s_sq/np.sum(s_sq)) < prop_var
    r_X = np.logical_and(dx,px)
    u = 

    return rates


