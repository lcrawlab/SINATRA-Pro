#!/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl
from matplotlib import cm

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'


data = []
for i in range(1,101):
    raw = np.loadtxt('ec_sphere/WT/WT_chimera_%d.dat'%i)
    x = raw[0]
    
    n_direction = raw.shape[0]-1
    n_length = x.size
    print(n_direction,n_length)

    fig, ax1 = plt.subplots(figsize=(8,6))
    colors = [cm.rainbow(x) for x in np.linspace(0,1,n_direction)]

    for idir in range(n_direction):
        y = raw[1+idir]
        ax1.plot(x,y,linestyle='-',color=colors[idir])

    ax1.set_xticks(np.arange(-30,32,10))
    xml = MultipleLocator(2)
    ax1.xaxis.set_minor_locator(xml)

    ax1.set_yticks(np.arange(-100,401,20))
    yml = MultipleLocator(5)
    ax1.yaxis.set_minor_locator(yml)

    ax1.set_xlabel(r'radius ($\mathrm{\AA}$)',fontsize=18)
    ax1.set_ylabel('EC',fontsize=18)
    ax1.tick_params(axis='y',labelsize=14)
    ax1.tick_params(axis='x',labelsize=14)

    ax1.set_xlim(-32,32)
    ax1.set_ylim(-20,400)
    #ax1.legend(loc='upper right',fontsize=16)

    plt.tight_layout()
    plt.savefig('ec_direction_%d.pdf'%i)
    #plt.show()
    plt.close()

