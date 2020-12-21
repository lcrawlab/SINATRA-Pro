#!/bin/python3
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'

directory='.'

mutant = "L67S"
pfix = "L67S"

ref = mda.Universe('L67S/serial/prot_init.gro')#%mutant)
traj = mda.Universe('%s/serial/prot_init.gro'%mutant,'%s/serial/prot_pbc_center.xtc'%mutant)

refCA = ref.select_atoms('name CA')
trajCA = traj.select_atoms('name CA')

rmsds = []
t = []

com_refCA = refCA.center_of_mass()
ref0 = refCA.positions - com_refCA
for ts in traj.trajectory:
    traj0 = trajCA.positions - trajCA.center_of_mass()
    R, rmsdval = align.rotation_matrix(traj0, ref0)
    trajCA.translate(-trajCA.center_of_mass())
    trajCA.rotate(R)
    trajCA.translate(com_refCA)
    t.append(ts.time)
    rmsds.append(rmsd(trajCA.positions,refCA.positions))

print(rmsds[0])    
print(np.mean(rmsds),np.std(rmsds))

np.savetxt(directory+'/rmsd_time_%s.dat'%pfix,np.transpose((t,rmsds)),fmt='%.6f')

#t,rmsds = np.loadtxt(directory+'/rmsd_time_WT.dat',unpack=True)

nbin = 161
bins = np.linspace(0,8,nbin)
hist, bin_edges = np.histogram(rmsds,bins=bins,density=True)
x = bins[:-1]+0.5*(bins[1]-bins[0])
y = hist
np.savetxt(directory+'/rmsd_hist_%s.dat'%pfix,np.transpose((x,y)),fmt='%.6f')

#t,rmsds = np.loadtxt(directory+'/rmsd_time.dat',unpack=True)

fig = plt.figure(figsize=(8,6))
ax1 = plt.subplot(1,1,1)

ax1.plot(t,rmsds,label='WT',color='black')

ax1.set_xticks(np.arange(0,100001,10000))
xml = MultipleLocator(2000)
ax1.xaxis.set_minor_locator(xml)
ax1.set_yticks(np.arange(0,10.0,1.0))
yml = MultipleLocator(0.5)
ax1.yaxis.set_minor_locator(yml)
ax1.set_xlabel('time (ps)',fontsize=18)
ax1.set_ylabel(r'RMSD ($\mathrm{\AA}$)',fontsize=18)
ax1.tick_params(axis='y',labelsize=14)
ax1.tick_params(axis='x',labelsize=14)

ax1.set_xlim(0,100000)
ax1.set_ylim(0,6.0)
#ax1.legend(loc='upper right',fontsize=16)

plt.tight_layout()
plt.savefig(directory+'/rmsd_time_%s.pdf'%pfix)
plt.show()
plt.close()

fig = plt.figure(figsize=(8,6))
ax1 = plt.subplot(1,1,1)

ax1.plot(bins[:-1]+0.5*(bins[1]-bins[0]),hist,color='black',marker='o')

ax1.set_xticks(np.arange(0.0,8.0,0.5))
xml = MultipleLocator(0.1)
ax1.xaxis.set_minor_locator(xml)
yml = MultipleLocator(0.1)
ax1.yaxis.set_minor_locator(yml)
ax1.set_xlabel(r'RMSD ($\mathrm{\AA}$)',fontsize=18)
ax1.set_ylabel('Normalized frequencies',fontsize=18)
ax1.tick_params(axis='y',labelsize=14)
ax1.tick_params(axis='x',labelsize=14)

ax1.set_xlim(0,8)
ax1.set_yticks(np.arange(0,4.1,0.2))
ax1.set_ylim(0,1)
#ax1.legend(loc='upper right',fontsize=16)

plt.tight_layout()
plt.savefig(directory+'/rmsd_dist_%s.pdf'%pfix)
plt.show()
plt.close()


