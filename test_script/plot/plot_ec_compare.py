#!/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'

datas = []
mutants = ['WT','WT chopped']
colors = ['blue','red']
ecfiles = ['sect_WT_65_213_4_15_01_4_25_python_normboth_mod_70_213.txt','sect_WT_65_213_mod_70_213_4_15_01_4_25_python_normboth_WT.txt']

length = 25

x = np.arange(-1.0,1.01,2./length)[1:]
for mutant, ecfile in zip(mutants,ecfiles):
    raw = np.loadtxt(ecfile)
    n_directions = int(raw.shape[1]/length)
    data = []
    for ifile in range(raw.shape[0]):
        row = []
        for i in range(n_directions):
            row.append(list(raw[ifile,length*i:length*(i+1)]))
        data.append(row)

    n_file = len(data)
    n_direction = len(data[0])
    n_length = len(data[0][0])
    print(n_file,n_direction,n_length)
    datas.append(data)

for idir in range(n_direction):
    fig, ax1 = plt.subplots(figsize=(8,6))
   
    for idata in range(len(datas)):
        data = datas[idata]
        color = colors[idata]
        mutant = mutants[idata]
        #inliers = np.loadtxt('inliers_%s.txt'%mutant).astype(int)
        ys = []
        for ifile in range(n_file): 
            y = data[ifile][idir]
            #y = np.append([0],y)
            ys.append(y)
            ax1.plot(x,y,linestyle='--',color=color,alpha=0.2)
        
        ys = np.array(ys)
        mean = np.mean(ys,axis=0)
        std = np.std(ys,axis=0)
        
        ax1.plot(x,mean,linestyle='-',color=color)
        ax1.fill_between(x,mean-std,mean+std,color=color,alpha=0.4)

    ax1.set_xticks(np.arange(-1.0,1.1,0.5))
    xml = MultipleLocator(0.1)
    ax1.xaxis.set_minor_locator(xml)

    #ax1.set_yticks(np.arange(0,15000,500))
    ax1.set_yticks(np.arange(-2,2,1))
    yml = MultipleLocator(0.1)
    ax1.yaxis.set_minor_locator(yml)

    ax1.set_xlabel('radius (unit sphere)',fontsize=18)
    ax1.set_ylabel('SEC',fontsize=18)
    ax1.tick_params(axis='y',labelsize=14)
    ax1.tick_params(axis='x',labelsize=14)

    ax1.set_xlim(-1.0,1.0)
    ax1.set_ylim(-1.5,1.5)
    
    #ax1.set_yticks(np.arange(-0.5,1.6,0.5))
    #yml = MultipleLocator(0.1)
    #ax1.set_ylim(-0.5,1.5)

    #ax1.legend(loc='upper right',fontsize=16)

    plt.tight_layout()
    plt.savefig('plot/sect_WT_65_213_mod_70_213_norm_25/sect_curve_%d_band.pdf'%idir)
    #plt.show()
    plt.close()

