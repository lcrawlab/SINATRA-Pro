#!/bin/python3

from mesh import *
from euler import *
from directions import *
from reconstruction import *

import matplotlib.pyplot as plt
import numpy as np
import os, sys

#from scipy.special import logsumexp
import MDAnalysis as mda

nf = 50
method = "dect"
protA = "WT_65_213"
protB = "R164S"

nframe = 101

meshProtein = mesh()

#directions = generate_equidistributed_cones(n_directions=15, directions_per_cone=4, cap_radius = 0.1, hemisphere=False)
#np.savetxt('directions_cone.txt',directions)
directions = np.loadtxt('directions_cone_15.txt')

#np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

#pip = np.loadtxt('pip_%s_%s_%s_%d_30.txt'%(protA,protB,method,nf))
#raw = np.loadtxt('pip_WT_R164S_65_213_dect_4_15_01_4_50.txt')
raw = np.loadtxt('rates_dect_WT_R164S_65_213_4_15_01_4_50_ESS_001.txt',usecols=1)
#raw = np.loadtxt('rates_dect_WT_R164S_65_213_50_ESS_001.txt',usecols=1)
not_vacuum = np.loadtxt('notvacuum_dect_WT_R164S_65_213_4_15_01_4_50_norm.txt')
pip = np.zeros(not_vacuum.size,dtype=float)
j = 0
for i in range(not_vacuum.size):
    if not_vacuum[i]:
        pip[i] = raw[j]
        j += 1

meshProtein.read_mesh_file(filename='mesh_noh/%s_4/WT_4_frame0.msh'%protA)
outProb = np.zeros((nframe,meshProtein.vertices.shape[0]),dtype=float)
for i in range(nframe):
    frame = 100*i
    sys.stdout.write('Reconstructing for Frame %d...\r'%frame)
    sys.stdout.flush()
    meshProtein.read_mesh_file(filename='mesh_noh/%s_4/WT_4_frame%d.msh'%(protA,frame))
    #outProb[i,:] = reconstruct_by_sorted_threshold(meshProtein, directions, pip, n_filtration = nf, n_direction_per_cone = 4, ball_radius = 1.0, by_rank = False)
    outProb[i,:] = reconstruct_by_sorted_threshold_new(meshProtein, directions, pip, n_filtration = nf, n_direction_per_cone = 4, ball_radius = 1.0, by_rank = False)
    #outProb[i,:] = reconstruct_by_average(meshProtein, directions, pip, n_filtration = nf, n_direction_per_cone = 4, ball_radius = 1.0)

averageProb = np.average(outProb,axis=0)
np.savetxt('prob_rate_sthrval_vert_%s_%s_%s_%d_15.txt'%(protA,protB,method,nf),averageProb)

prob = np.loadtxt('prob_rate_sthrval_vert_%s_%s_%s_%d_15.txt'%(protA,protB,method,nf))
u = mda.Universe('pdb_noh/%s/WT_frame5000.pdb'%protA)
protein = u.select_atoms('protein and not type H')
u.add_TopologyAttr('tempfactors')
print(np.amin(prob),np.amax(prob))
y = prob*100
#y = prob*2500
#y = prob/np.amax(prob)*25
print(np.amin(y),np.amax(y))
protein.tempfactors = y
protein.write("%s_%s_%s_%d_15_rate_sthrval.pdb"%(protA,protB,method,nf))

