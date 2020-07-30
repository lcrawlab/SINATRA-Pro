#!/bin/python3

from mesh import *
from euler import *
from directions import *

import matplotlib.pyplot as plt
import numpy as np
import os

meshProtein = mesh()

directions = generate_equidistributed_cones(n_directions=10, directions_per_cone=10, cap_radius = 0.3, hemisphere=False)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for direction in directions:
    ax.scatter(direction[0],direction[1],direction[2])
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)
plt.show()

np.savetxt('directions_cone.txt',directions)
#directions = np.loadtxt('directions_cone.txt')

for i in range(1,101):
    meshProtein.read_obj_file(filename='rssb/mesh/WT/WT_chimera_%d.obj'%i)#,read_edges_from_file=True,edges_filename='edges/WT/WT_chimera_%d.txt'%i)
    #meshProtein.save_edges_file('edges/WT/WT_chimera_%d.txt'%i)
    outfilename = 'rssb/ec_sphere/WT/WT_chimera_%d.dat'%i
    if os.path.exists(outfilename):
        open(outfilename, 'w').close()
    with open(outfilename,'ab') as f:
        for j, direction in enumerate(directions):
            vf = compute_vertex_function(meshProtein,direction)
            t, ec = compute_ec_curve(meshProtein,vf,curve_length=100,ball_radius=50.0,standardized=True,ec_type="EC")
            if j == 0:
                np.savetxt(f,[t],fmt="%.8f",delimiter=" ")
            np.savetxt(f,[ec],fmt="%.8e",delimiter=" ")

    


