#!/bin/python3

#from mesh import *
#from euler import *
#from directions import *
from gp import *

#import matplotlib.pyplot as plt
import numpy as np
import os, sys

x = np.loadtxt("dect_WT_R164S_65_213_4_15_01_4_25_norm.txt")
N = x.shape[0]
y = np.zeros(N,dtype=int)
y[:int(N/2)].fill(-1)
y[int(N/2):].fill(1)

rates = find_rate_variables_with_other_sampling_methods(x,y,bandwidth=0.01,sampling_method="ESS")

