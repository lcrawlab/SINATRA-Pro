#!/bin/python3

from mesh import *
from euler import *
from directions import *

import matplotlib.pyplot as plt
import numpy as np
import os, sys

mutants = ["WT","R164S","R164G"]

meshProtein = mesh()

#directions = generate_equidistributed_cones(n_directions=30, directions_per_cone=4, cap_radius = 0.1, hemisphere=False)
#np.savetxt('directions_cone_30.txt',directions)
directions = np.loadtxt('directions_cone.txt')

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

for mutant in mutants:
    ecs = []
    for i in range(101):
        frame = i*100
        sys.stdout.write('Calculating EC for %s Frame %d...\r'%(mutant,frame))
        sys.stdout.flush()
        meshProtein.read_mesh_file(filename='mesh_noh/%s_65_213_4/%s_4_frame%d.msh'%(mutant,mutant,frame))
        t, ec = compute_ec_curve(meshProtein,directions,n_filtration=100,ball_radius=1.0,ec_type="SECT",include_faces=False)
        ecs.append(ec.flatten())
    ecs = np.array(ecs)
    np.savetxt('sect_nofaces_%s_65_213_4_15_01_4_100.txt'%mutant,ecs,fmt='%.3f')

