#!/bin/python3

from mesh import *
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import to_rgb
from scipy.spatial import Delaunay

from numba import jit

@jit(nopython=True)
def calc_density_grid(density,dx,dy,dz,n_grid,alpha=1.0,lmbda=1.0):
    for i in range(n_grid):
        for j in range(n_grid):
            for k in range(n_grid):
                density[i,j,k] += alpha/(dx[i]*dx[i]+dy[j]*dy[j]+dz[k]*dz[k]+lmbda)
    return

def calc_density(meshA,radius=1.0,n_grid=11,alpha=1.0,lmbda=1.0):
    p = np.linspace(-radius,radius,n_grid)
    density = np.zeros((n_grid,n_grid,n_grid),dtype=float) 
    for i in range(meshA.n_vertex):
        x = meshA.vertices[i]
        dx = x[0]-p
        dy = x[1]-p
        dz = x[2]-p
        calc_density_grid(density,dx,dy,dz,n_grid)
    return density

'''
def calc_density_new(meshA,radius=1.0,n_grid=11,alpha=1.0,lmbda=1.0):
    p = np.linspace(-radius,radius,n_grid)
    density = np.zeros((n_grid,n_grid,n_grid),dtype=float) 
    for i in range(meshA.n_vertex):
        x = meshA.vertices[i]
        dx = x[0]-p
        dy = x[1]-p
        dz = x[2]-p
        dx = dx*dx
        dy = dy*dy
        dz = dz*dz
        dx = np.transpose(np.tile(np.tile(dx,n_grid),n_grid).reshape(n_grid,n_grid,n_grid),(2,1,0))
        dy = np.transpose(np.tile(np.tile(dy,n_grid),n_grid).reshape(n_grid,n_grid,n_grid),(1,2,0))
        dz = np.transpose(np.tile(np.tile(dz,n_grid),n_grid).reshape(n_grid,n_grid,n_grid),(0,1,2))
        density += alpha/(dx+dy+dz+lmbda)
    return density
'''

def calc_level(density,n_level=10):
    ymax = np.amax(density)
    bins = np.linspace(0,ymax,n_level+1)
    result = np.digitize(density,bins)
    return result

def triangulation(radius,n_grid,levelset,n_level,cutoff=1,calc_heat=False,seperate_off_file=True,filename="output"):
    if calc_heat:
        heat = []
        colors = cm.rainbow(np.linspace(0,1,n_level))
        spectrum = np.floor(np.array([to_rgb(a) for a in colors])*255).astype(int)
    x = np.linspace(-radius,radius,n_grid)
    grid = np.array(np.meshgrid(x,x,x))
    meshB = mesh()
    for i in range(cutoff,n_level):
        if seperate_off_file:
            meshB = mesh()
            if calc_heat:
                heat = []
        level = levelset == i
        points = np.transpose(grid[:,level])
        if len(meshB.vertices) == 0:
            meshB.vertices = points
        else:
            meshB.vertices = np.vstack((meshB.vertices,points))
        if len(points) > 0:
            tri = Delaunay(points)
            if len(meshB.faces) == 0:
                meshB.faces = tri.simplices
            else:
                meshB.faces = np.vstack((meshB.faces,tri.simplices))
            if calc_heat:
                for j in range(tri.simplices.shape[0]):
                    heat.append(spectrum[i])
        if seperate_off_file:
            if calc_heat:
                heat = np.array(heat)
                meshB.write_off_file_heat(heat=heat,filename=filename+"_%d_heat.off"%i)
            else:
                meshB.write_off_file(filename=filename+"_%d.off"%i)
    if calc_heat:
        return np.array(heat)
    return


