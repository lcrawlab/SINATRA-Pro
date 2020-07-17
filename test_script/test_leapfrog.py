#!/bin/python3


import numpy as np
from matplotlib import pyplot as plt

init_x = np.array([0.0,0.0])
init_v = np.array([10.0,10.0])

def a(x):
    return np.array([0,-10])

plt.figure()

dt = 0.01
v = init_v + a(init_x)*dt*.5
x = init_x + v*dt

for i in range(100):
    v = v + a(x)*dt
    x = x + v*dt
    plt.scatter(x[0],x[1])
plt.show()




