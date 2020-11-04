#!/bin/python3

from mesh import *
from euler import *
from directions import *

import matplotlib.pyplot as plt
import numpy as np
import os, sys

mutants = ["WT"]#,"R164S","R164G"]

level = 9
n_filtration = 25
n_direction = 15
method = "ECT"

directions = generate_equidistributed_cones(n_directions=n_direction, directions_per_cone=4, cap_radius = 0.1, hemisphere=False)
np.savetxt('directions_cone_%d.txt'%n_direction,directions)
directions = np.loadtxt('directions_cone_%d.txt'%n_direction)

for mutant in mutants:
    ecs = []
    for i in range(101):
        frame = i*100
        sys.stdout.write('Calculating EC for %s Frame %d...\r'%(mutant,frame))
        sys.stdout.flush()
        meshProtein = mesh()
        #meshProtein.read_off_file(filename='mesh/%s_%d/%s_frame%d.off'%(mutant,level,mutant,frame))
        #meshProtein.write_mesh_file(edges=meshProtein.edges,faces=meshProtein.faces,filename='mesh/%s_%d/%s_frame%d.msh'%(mutant,level,mutant,frame))
        meshProtein.read_mesh_file(filename='mesh/%s_%d/%s_frame%d.msh'%(mutant,level,mutant,frame))
        t, ec = compute_ec_curve(meshProtein,directions,n_filtration=n_filtration,ball_radius=1.0,ec_type=method,include_faces=True)
        ecs.append(ec.flatten()) 
    ecs = np.array(ecs)
    np.savetxt('%s_%s_65_213_%d_%d_01_4_%d.txt'%(method,mutant,level,n_direction,n_filtration),ecs,fmt='%d')


