#!/bin/python3

from mesh import *
from euler import *
from directions import *

import matplotlib.pyplot as plt
import numpy as np
import os
import time


#directions = generate_equidistributed_cones(n_directions=100, directions_per_cone=1, hemisphere=False)
#np.savetxt('directions_1.txt',directions)
directions = np.loadtxt('directions_1.txt')

#os.mkdir('ec_sphere_3/')
#os.mkdir('ec_sphere_3/R164G/')

for i in range(1,101):
    meshProtein = mesh()
    meshProtein.read_obj_file(filename='mesh/R164G/R164G_chimera_%d.obj'%i)#,read_edges_from_file=True,edges_filename='edges/R164G/R164G_chimera_%d.txt'%i)
    #meshProtein.save_edges_file('edges/R164G/R164G_chimera_%d.txt'%i)
    outfilename = 'ec_sphere_3/R164G/R164G_chimera_%d.dat'%i
    if os.path.exists(outfilename):
        open(outfilename, 'w').close()

    directions = directions[:20]

    start = time.time()

    with open(outfilename,'ab') as f:
        t, ecs = compute_ec_curve(meshProtein,directions,n_filtration=100,ball_radius=40.0,standardized=True,ec_type="EC")
        end = time.time()
        print(end-start)


        np.savetxt(f,[t],fmt="%.8f",delimiter=" ")
        for ec in ecs:
            np.savetxt(f,[ec],fmt="%.8e",delimiter=" ")

        exit()


