#!/bin/python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'

directory='.'

fig, ax1 = plt.subplots(figsize=(8,6))
#ax1 = plt.subplot(1,1,1)

ec_directory = 'ec_curve/'
for i in range(1,21):
    data = np.loadtxt(ec_directory+'SECT_%d.dat'%i)
    x = data[0]
    y = data[1]
    ax1.plot(x,y)#,label='WT',color='black')

ax1.set_xticks(np.arange(-50,51,10))
xml = MultipleLocator(2)
ax1.xaxis.set_minor_locator(xml)

ax1.set_yticks(np.arange(-500,4001,500))
yml = MultipleLocator(100)
ax1.yaxis.set_minor_locator(yml)

ax1.set_xlabel(r'radius ($\mathrm{\AA}$)',fontsize=18)
ax1.set_ylabel('SECT',fontsize=18)
ax1.tick_params(axis='y',labelsize=14)
ax1.tick_params(axis='x',labelsize=14)

ax1.set_xlim(-50,50)
ax1.set_ylim(-500,4000)
#ax1.legend(loc='upper right',fontsize=16)

plt.tight_layout()
plt.savefig(directory+'/sect_curve.pdf')
plt.show()
plt.close()

