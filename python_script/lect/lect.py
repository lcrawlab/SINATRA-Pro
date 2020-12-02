#!/bin/python3

from mesh import *
import numpy as np
from numba import jit
from matplotlib import pyplot as plt
from matplotlib import cm
from matplotlib.colors import to_rgb
from scipy.spatial import Delaunay
from scipy.spatial.distance import pdist

@jit(nopython=True)
def calc_density_grid(density,dx,dy,dz,n_grid,alpha=1.0,lmbda=1.0):
    for i in range(n_grid):
        for j in range(n_grid):
            for k in range(n_grid):
                density[i,j,k] += alpha/(dx[i]*dx[i]+dy[j]*dy[j]+dz[k]*dz[k]+lmbda)
    return

## calc_density_restriction allows a boolean matrix input that specific which grid to calculate
## e.g. grid point in the sphere = TRUE, grid point outside the sphere = FALSE
@jit(nopython=True)
def calc_density_restriction(density,calc_this,dx,dy,dz,n_grid,alpha=1.0,lmbda=1.0):
    for i in range(n_grid):
        for j in range(n_grid):
            for k in range(n_grid):
                if calc_this[i,j,k]:
                    density[i,j,k] += alpha/(dx[i]*dx[i]+dy[j]*dy[j]+dz[k]*dz[k]+lmbda)
    return

## Density calculation by summing contribution from all atoms in the protein
def calc_density(meshA,radius=1.0,n_grid=11,alpha=1.0,lmbda=1.0,calc_this=[],restriction=False):
    p = np.linspace(-radius,radius,n_grid)
    density = np.zeros((n_grid,n_grid,n_grid),dtype=float)
    if restriction:
        for i in range(meshA.n_vertex):
            x = meshA.vertices[i]
            dx = x[0]-p
            dy = x[1]-p
            dz = x[2]-p
            calc_density_restriction(density,calc_this,dx,dy,dz,n_grid,alpha,lmbda)
    else:
        for i in range(meshA.n_vertex):
            x = meshA.vertices[i]
            dx = x[0]-p
            dy = x[1]-p
            dz = x[2]-p
            calc_density_grid(density,dx,dy,dz,n_grid,alpha,lmbda)
    return density

def calc_level(density,n_level=10):
    ymax = np.amax(density)
    bins = np.linspace(0,ymax,n_level+1)
    result = np.digitize(density,bins)
    return result

def triangulation(radius,n_grid,levelset,n_level,cutoff=1,n_shapecutoff=5,calc_heat=False,seperate_off_file=True,mutant="WT",frame=0):
    if calc_heat:
        heat = []
        colors = cm.rainbow(np.linspace(0,1,n_level))
        spectrum = np.floor(np.array([to_rgb(a) for a in colors])*255).astype(int)
    x = np.linspace(-radius,radius,n_grid)
    gridsize = x[1]-x[0]
    shapecutoff = gridsize*n_shapecutoff
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
            for face in tri.simplices:
                vertices = points[face]
                dist = pdist(vertices)
                if np.any(dist > shapecutoff):
                    continue
                else:
                    meshB.faces.append(face)
                    if calc_heat:
                        heat.append(spectrum[i])
        if seperate_off_file:
            meshB.faces = np.array(meshB.faces)
            if calc_heat:
                heat = np.array(heat)
                #meshB.write_off_file_heat(heat=heat,filename="mesh/%s_%d/%s_frame%d.off"%(mutant,i,mutant,frame))
            else:
                meshB.write_off_file(filename="mesh/%s_%d/%s_frame%d.off"%(mutant,i,mutant,frame))
                #meshB.write_off_file(filename="mesh/%s_%d/%s_%d_frame%d.off"%(mutant,i,mutant,n_shapecutoff,frame))
    if not seperate_off_file:
        meshB.faces = np.array(meshB.faces)
    if calc_heat:
        return np.array(heat)
    return

