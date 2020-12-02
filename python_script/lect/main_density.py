#!/bin/python3

from mesh import *
from lect import *

import matplotlib.pyplot as plt
import numpy as np
import os, sys
import MDAnalysis as mda

mutants = ["WT","R164S"]
n_grid = 21
n_sample = 101

###
# Calculate maxium radius among all structures
###

radii = []
for mutant in mutants:
    for i in range(101):
        frame = i*100
        trajwhole = mda.Universe('pdb_noh/%s_65_213/%s_frame%d.pdb'%(mutant,mutant,frame)).select_atoms('protein and not type H')
        sys.stdout.write('Calculating grid size for %s for Frame %d...\r'%(mutant,frame))
        sys.stdout.flush()
        meshA = mesh()
        meshA.vertices = trajwhole.positions
        meshA.n_vertex = meshA.vertices.shape[0]
        radius = np.amax(np.linalg.norm(meshA.vertices,axis=1))
        radii.append(radius)
radius = np.amax(radii)
#np.savetxt("max_radius.dat",[radius],fmt="%.6f")
#radius = np.loadtxt("max_radius.dat",dtype=float)
#print(radius)

###
# Generate a spherical grid with the maximum radius
###

gridsize = radius / (n_grid-1)
calc_this = np.zeros((n_grid,n_grid,n_grid),dtype=bool)
middle = int((n_grid-1)/2)+1
rsq = ((n_grid-1)/2)**2
position = []
n_point = 0
p = np.linspace(-radius,radius,n_grid)
indices = np.zeros((n_grid,n_grid,n_grid),dtype=int)
for i in range(n_grid):
    for j in range(n_grid):
        for k in range(n_grid):
            r = (i-middle)**2+(j-middle)**2+(k-middle)**2
            if r <= rsq:
                calc_this[i,j,k] = True
                indices[i,j,k] = n_point
                position.append([p[i],p[j],p[k]])
                n_point += 1
np.save('indices_%d.npy'%n_grid,indices)
np.savetxt('position_%d.dat'%n_grid,position,fmt="%.6f")


###
# Calculate the density field at each grid point for each structure
###
for mutant in mutants:
    density_matrix = np.zeros((n_sample,n_point),dtype=float)
    for i in range(n_sample):
        frame = i*100
        trajwhole = mda.Universe('pdb_noh/%s_65_213/%s_frame%d.pdb'%(mutant,mutant,frame)).select_atoms('protein and not type H')
        sys.stdout.write('Calculating grid size for %s for Frame %d...\r'%(mutant,frame))
        sys.stdout.flush()
        meshA = mesh()
        meshA.vertices = trajwhole.positions
        meshA.n_vertex = meshA.vertices.shape[0]
        density = calc_density(meshA=meshA,radius=radius,n_grid=n_grid,alpha=1.0,lmbda=1.0,calc_this=calc_this,restriction=True)
        density = density[calc_this]
        density = density.flatten()
        density_matrix[i] = density
    np.savetxt('density_%d_%s.dat'%(n_grid,mutant),density_matrix,fmt="%.6f")

