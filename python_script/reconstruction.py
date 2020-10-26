#!/bin/python3

import numpy as np
from scipy.special import logsumexp

# Reconstruction

def reconstruct_on_vertices(mesh, directions, pip, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0, standardized = True):
    outProb = np.zeros(mesh.vertices.shape[0],dtype=float)
    n_direction = directions.shape[0]
    n_cone = int(n_direction / n_direction_per_cone)
    for i in range(n_cone):
        coneProb = np.zeros(mesh.vertices.shape[0],dtype=float)
        for j in range(n_direction_per_cone):
            k = i*n_direction_per_cone+j
            vertex_function = np.dot(mesh.vertices,directions[k])
            if standardized:
                radius = np.linspace(-ball_radius,ball_radius,n_filtration)
            else:
                radius = np.linspace(np.amin(vertex_function),np.amax(vertex_function),n_filtration)        
            filtration = np.digitize(vertex_function,radius)
#            coneProb += pip[k*n_filtration+filtration]
#        coneProb /= n_direction_per_cone
#        outProb += np.log(coneProb)
#    outProb -= logsumexp(outProb)
            outProb += pip[k*n_filtration+filtration]
    outProb /= n_direction
    #outProb /= np.sum(outProb)
    return outProb      

def reconstruct_on_vertices_multiply(mesh, directions, pip, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0, standardized = True):
    outProb = np.zeros(mesh.vertices.shape[0],dtype=float)
    n_direction = directions.shape[0]
    n_cone = int(n_direction / n_direction_per_cone)
    for i in range(n_cone):
        coneProb = np.zeros(mesh.vertices.shape[0],dtype=float)
        for j in range(n_direction_per_cone):
            k = i*n_direction_per_cone+j
            vertex_function = np.dot(mesh.vertices,directions[k])
            if standardized:
                radius = np.linspace(-ball_radius,ball_radius,n_filtration)
            else:
                radius = np.linspace(np.amin(vertex_function),np.amax(vertex_function),n_filtration)        
            filtration = np.digitize(vertex_function,radius)
            coneProb += pip[k*n_filtration+filtration]
        coneProb /= n_direction_per_cone
        outProb += np.log(coneProb)
    outProb -= logsumexp(outProb)
    return outProb      


        
