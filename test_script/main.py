#!/bin/python3

from mesh import *
from euler import *
from directions import *

import matplotlib.pyplot as plt
import numpy as np

## Test initialization and import topology files

#mesh = mesh()
#mesh.read_off_file(filename='data/bean_off_file.off',read_edges_from_file=true,edges_filename='edges.txt')

#mesh.read_obj_file(filename='data/1xpb_chimera.obj',read_edges_from_file=true,edges_filename='1xpb_edges.txt')
#mesh.save_edges_file('edges_1xpb.txt')
#mesh.centering()

#print("# Vertices = %d"%mesh.n_vertex)
#print("# Edges = %d"%mesh.n_edge)
#print("# Faces = %d"%mesh.n_face)

'''
## Testing opposite direction of filtration
fig, ax = plt.subplots()
vf = compute_vertex_function(mesh,[-1,1,0])
t, ec = compute_ec_curve(mesh,vf,curve_length=1000,ball_radius=40.0,standardized=True,ec_type="EC")
ax.plot(t,ec,'b-')
#vf = compute_vertex_function(mesh,[1,-1,0])
#t, ec = compute_ec_curve(mesh,vf,curve_length=1000,ball_radius=40.0,standardized=True,ec_type="EC")
#ax.plot(t,ec,'r-')
#t, ec = compute_ec_curve_test(mesh,vf,curve_length=1000,ball_radius=40.0,standardized=True,ec_type="EC")
#ax.plot(t,ec,'r--')
plt.show()
'''

## Test generated cone directions, cone script not completely checked
#directions = generate_equidistributed_cones(n_directions=20, directions_per_cone=1, hemisphere=False)

'''
## Visualize directions over sphere
fig = plt.figure(figsize=plt.figaspect(0.5)*1.5)
ax = fig.gca(projection='3d')
for direction in directions:
    ax.scatter(direction[0],direction[1],direction[2])
plt.show()
'''
mesh = mesh()

ec_directory = 'ec_curve/WT/'
direction = [1,0,0]

fig, ax = plt.subplots()
for i in range(1,6):
    mesh.read_obj_file(filename='data/WT/WT_%d.obj'%i)#,read_edges_from_file=True,edges_filename='1xpb_edges.txt')
    #mesh.save_edges_file('edges_1xpb.txt')
    mesh.centering()

    vf = compute_vertex_function(mesh,direction)
    t, ec = compute_ec_curve(mesh,vf,curve_length=100,ball_radius=50.0,standardized=True,ec_type="SECT")
    np.savetxt(ec_directory+'SECT_WT_%d.dat'%i,(t,ec))
    ax.plot(t,ec,'-')

plt.savefig('SECT_WT.pdf')
plt.show()


