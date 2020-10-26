#!/bin/python3

import numpy as np
from matplotlib import pyplot as plt

import csv

heat = []
with open("remote_data/vert_heat_control.txt","r") as f:
    p = csv.reader(f)
    for row in p:
        heat.append(row[0])

vertices = []
with open("remote_data/complex_vertices.txt","r") as f:
    p = csv.reader(f,delimiter=" ")
    for row in p:
        vertices.append(row[1:])

vertices = np.array(vertices[1:]).astype(float)

print(len(heat))
print(vertices.shape)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.scatter(vertices[:,0],vertices[:,1],vertices[:,2],c=heat,s=1)
plt.savefig('heat_mesh_control.jpg')
plt.show()

