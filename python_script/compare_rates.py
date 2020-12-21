#!/bin/python3

from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import linregress 

x = np.loadtxt("data/rates_R.txt",usecols=1)
#y = np.loadtxt("data/rates_python_Ress_Rsvd.txt")
y = np.loadtxt("data/rates_python.txt")

#plt.plot(x)
#plt.plot(y)
#plt.show()

slope, intercept, r_value, p_value, std_err = linregress(x, y)
plt.plot(x, intercept + slope*x, 'r')
plt.scatter(x,y)
plt.xlabel("RATEs (R)")
plt.ylabel("RATEs (python)")
plt.show()


