#!/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl
from matplotlib import cm

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'

n_file = 100

ec = []
for i in range(n_file):
    raw = np.loadtxt('ec_sphere/WT/WT_chimera_%d.dat'%(i+1))
    ec.append(raw[1:,-1])

fig, ax1 = plt.subplots(figsize=(8,6))
for i in range(n_file):
    n_direction = len(ec[i])
    colors = [cm.rainbow(x) for x in np.linspace(0,1,n_direction)]
    for idir in range(n_direction):
        ax1.scatter([i]*n_direction,ec[i],color=colors)

ax1.set_xticks(np.arange(0,100,10))
xml = MultipleLocator(1)
ax1.xaxis.set_minor_locator(xml)

ax1.set_yticks(np.arange(0,1001,200))
yml = MultipleLocator(50)
ax1.yaxis.set_minor_locator(yml)

ax1.set_xlabel('frame #',fontsize=18)
ax1.set_ylabel('EC',fontsize=18)
ax1.tick_params(axis='y',labelsize=14)
ax1.tick_params(axis='x',labelsize=14)

ax1.set_xlim(0,100)
ax1.set_ylim(0,1000)
#ax1.legend(loc='upper right',fontsize=16)

plt.tight_layout()
plt.savefig('ec_total.pdf')
plt.close()

