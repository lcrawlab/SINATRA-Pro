#!/bin/python3

from mesh import *
from euler_parallel import *
from directions import *

import numpy as np
import os, sys

mutants = ["WT","R164S"]#,"R164G"]

meshProtein = mesh()

n_cone = 15
n_fil = 100

directions = generate_equidistributed_cones(n_directions=n_cone, directions_per_cone=4, cap_radius = 0.1, hemisphere=False)
np.savetxt('directions_cone_%d.txt'%n_cone,directions)
directions = np.loadtxt('directions_cone_%d.txt'%n_cone)

np.set_printoptions(formatter={'float': '{: 0.3f}'.format})

for mutant in mutants:
    ecs = []
    for i in range(101):
        frame = i*100
        sys.stdout.write('Calculating EC for %s Frame %d...\r'%(mutant,frame))
        sys.stdout.flush()
        meshProtein.read_mesh_file(filename='mesh_noh/%s_65_213_4/%s_4_frame%d.msh'%(mutant,mutant,frame))
        t, ec = compute_ec_curve_parallel(meshProtein,directions,n_filtration=n_fil,ball_radius=1.0,ec_type="DECT",include_faces=True,n_core=6)
        ecs.append(ec.flatten())
    ecs = np.array(ecs)
    np.savetxt('dect_%s_65_213_4_%d_01_4_%d.txt'%(mutant,n_cone,n_fil),ecs,fmt='%.3f')

