#!/bin/python3

from mesh import *
from euler import *
from directions import *
from reconstruction import *

import matplotlib.pyplot as plt
import numpy as np
import os, sys

from scipy.special import logsumexp

meshProtein = mesh()

#directions = generate_equidistributed_cones(n_directions=15, directions_per_cone=4, cap_radius = 0.1, hemisphere=False)
#np.savetxt('directions_cone.txt',directions)
directions = np.loadtxt('directions_cone.txt')

#np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

nf = 50

pip = np.loadtxt('pip_WT_65_213_R164S_ect_nofaces_%d.txt'%nf)
#pip[pip == 1] = np.amax(pip[pip < 1])
pip /= np.sum(pip)

nframe = 101

meshProtein.read_off_file(filename='off_noh/WT_65_213_4/WT_4_frame0.off')
outProb = np.zeros((nframe,meshProtein.vertices.shape[0]),dtype=float)
for i in range(nframe):
    frame = 100*i
    sys.stdout.write('Reconstructing for Frame %d...\r'%frame)
    sys.stdout.flush()
    meshProtein.read_off_file(filename='off_noh/WT_65_213_4/WT_4_frame%d.off'%frame)
    outProb[i,:] = reconstruct_on_vertices(meshProtein, directions, pip, n_filtration = nf, n_direction_per_cone = 4, ball_radius = 1.0, standardized = True)

#averageProb = logsumexp(outProb,axis=0)
averageProb = np.average(outProb,axis=0)
np.savetxt('prob_vert_WT_65_213_R164S_ect_nofaces_%d.txt'%nf,averageProb)

