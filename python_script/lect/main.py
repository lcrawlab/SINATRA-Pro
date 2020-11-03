#!/bin/python3

from mesh import *
from lect import *

import matplotlib.pyplot as plt
import numpy as np
import os, sys
import MDAnalysis as mda

mutant = "WT"


np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

frame = 0

trajwhole = mda.Universe('../pdb_noh/%s_65_213/%s_frame%d.pdb'%(mutant,mutant,frame)).select_atoms('protein and not type H')
#    sys.stdout.write('Constructing topology for %s for Frame %d...\r'%(mutant,frame))
#    sys.stdout.flush()

n_grid = 51
n_level = 10

meshA = mesh()
meshA.vertices = trajwhole.positions
meshA.n_vertex = meshA.vertices.shape[0]
radius = np.amax(np.linalg.norm(meshA.vertices,axis=1))
density = calc_density(meshA=meshA,radius=radius,n_grid=n_grid)
levelset = calc_level(density,n_level=n_level)

meshB = mesh()
heat = triangulation(meshB,radius,n_grid,levelset,n_level=n_level,cutoff=5,calc_heat=True)
heat = np.array(heat)
meshB.write_off_file_heat(heat,filename='output_5.off')


