#!/bin/python3

from mesh import *
from euler import *
from directions import *

import matplotlib.pyplot as plt
import numpy as np

## Test initialization and import OFF files
mesh = mesh()
#meshCube.read_off_file(filename='data/1xpb_pymol.off')

#mesh.read_off_file(filename='data/bean_off_file.off',read_edges_from_file=True,edges_filename='edges.txt')
mesh.read_obj_file(filename='data/1xpb_chimera.obj',read_edges_from_file=True,edges_filename='1xpb_edges.txt')
#mesh.save_edges_file('edges_1xpb.txt')

print("# Vertices = %d"%mesh.n_vertex)
print("# Edges = %d"%mesh.n_edge)
print("# Faces = %d"%mesh.n_face)

ecs = []


    vf = compute_vertex_function(mesh,[1,0,0])
    t, ec = compute_ec_curve(mesh,vf,curve_length=100,ball_radius=20.0,standardized=False,ec_type="SECT")


fig, ax = plt.subplots()
ax.plot(t,ec,'-')
plt.show()


## Test generated cone directions

#points = generate_equidistributed_cones(n_directions=10, cap_radius=2, directions_per_cone=10)

##from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import

#fig = plt.figure(figsize=plt.figaspect(1))
#ax = fig.gca(projection='3d')

#for point in points:
#    ax.scatter(point[0],point[1],point[2])

#plt.show()


