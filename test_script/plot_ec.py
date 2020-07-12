#!/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'

data = []
for i in range(1,101):
    raw = np.loadtxt('ec_sphere/R164G/R164G_chimera_%d.dat'%i)
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
        if y[-1] < 300 and y[-1] > 100:
            ys.append(y)
        ax1.plot(x,y,linestyle='--',color='red',alpha=0.2)
    
    ys = np.array(ys)
    mean = np.mean(ys,axis=0)
    std = np.std(ys,axis=0)
    
    ax1.plot(x,mean,linestyle='-',color='red')
    ax1.fill_between(x,mean-std,mean+std,color='red',alpha=0.4)

    ax1.set_xticks(np.arange(-40,41,10))
    xml = MultipleLocator(2)
    ax1.xaxis.set_minor_locator(xml)

    ax1.set_yticks(np.arange(-200,1501,100))
    yml = MultipleLocator(20)
    ax1.yaxis.set_minor_locator(yml)

    ax1.set_xlabel(r'radius ($\mathrm{\AA}$)',fontsize=18)
    ax1.set_ylabel('EC',fontsize=18)
    ax1.tick_params(axis='y',labelsize=14)
    ax1.tick_params(axis='x',labelsize=14)

    ax1.set_xlim(-40,40)
    ax1.set_ylim(-50,200)
    #ax1.legend(loc='upper right',fontsize=16)

    plt.tight_layout()
    plt.savefig('plot/R164G/ec/ec_curve_%d.pdf'%idir)
    #plt.show()
    plt.close()

