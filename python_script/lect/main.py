#!/bin/python3

from mesh import *
from lect import *

import matplotlib.pyplot as plt
import numpy as np
import os, sys
import MDAnalysis as mda

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

trajwhole = mda.Universe('WT_frame0.pdb').select_atoms('protein and not type H')

n_grid = 51
n_level = 10

meshA = mesh()
meshA.vertices = trajwhole.positions
meshA.n_vertex = meshA.vertices.shape[0]
radius = np.amax(np.linalg.norm(meshA.vertices,axis=1))
density = calc_density(meshA=meshA,radius=radius,n_grid=n_grid)
levelset = calc_level(density,n_level=n_level)

meshB = mesh()
heat = triangulation(radius=radius,n_grid=n_grid,levelset=levelset,n_level=n_level,cutoff=1,calc_heat=True,seperate_off_file=True)
#meshB.write_off_file_heat(heat,filename='output_5.off')


