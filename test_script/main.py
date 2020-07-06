#!/bin/python3

from mesh import *
from euler import *
from directions import *

## Test initialization and import OFF files
meshCube = mesh()
meshCube.read_off_file(filename='data/bean_off_file.off')
#print(meshCube.vertices.shape)
#print(meshCube.edges.shape)
#print(len(meshCube.faces))

#vf = compute_vertex_function(meshCube,[1,0,0])
#output = compute_ec_curve(meshCube,vf,standardized=True,ec_type="EC")
#print(output)


## Test generated cone directions

#points = generate_equidistributed_cones(n_directions=10, cap_radius=2, directions_per_cone=10)

##from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
#import matplotlib.pyplot as plt
#import numpy as np

#fig = plt.figure(figsize=plt.figaspect(1))
#ax = fig.gca(projection='3d')

#for point in points:
#    ax.scatter(point[0],point[1],point[2])

#plt.show()


