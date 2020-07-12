#!/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'

data = []
for i in range(1,101):
    raw = np.loadtxt('ec_sphere/WT/WT_chimera_%d.dat'%i)
    if i == 1:
        x = raw[0]
    y = raw[1:]
    data.append(y)

n_file = len(data)
n_direction = len(data[0])
n_length = len(data[0][0])
print(n_file,n_direction,n_length)

for idir in range(n_direction):
    fig, ax1 = plt.subplots(figsize=(8,6))
    
    ys = []
    for ifile in range(n_file): 
        y = data[ifile][idir]
        ys.append(y)
        ax1.plot(x,y,linestyle='--',color='blue',alpha=0.2)
    
    ys = np.array(ys)
    mean = np.mean(ys,axis=0)
    std = np.std(ys,axis=0)
    
    ax1.plot(x,mean,linestyle='-',color='blue')
    ax1.fill_between(x,mean-std,mean+std,color='blue',alpha=0.4)

    ax1.set_xticks(np.arange(-30,32,10))
    xml = MultipleLocator(2)
    ax1.xaxis.set_minor_locator(xml)

    ax1.set_yticks(np.arange(-100,1501,100))
    yml = MultipleLocator(20)
    ax1.yaxis.set_minor_locator(yml)

    ax1.set_xlabel(r'radius ($\mathrm{\AA}$)',fontsize=18)
    ax1.set_ylabel('SECT',fontsize=18)
    ax1.tick_params(axis='y',labelsize=14)
    ax1.tick_params(axis='x',labelsize=14)

    ax1.set_xlim(-32,32)
    ax1.set_ylim(-100,1500)
    #ax1.legend(loc='upper right',fontsize=16)

    plt.tight_layout()
    #plt.savefig('ec_curve_%d.pdf'%i)
    plt.show()
    plt.close()
    break

