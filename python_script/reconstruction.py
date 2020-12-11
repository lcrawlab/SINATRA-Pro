#!/bin/python3

from mesh import *
import numpy as np
from scipy.special import logsumexp

# Reconstruction

def reconstruct_by_average(meshA, directions, pip, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0, standardized = True):
    outProb = np.zeros(meshA.vertices.shape[0],dtype=float)
    n_direction = directions.shape[0]
    n_cone = int(n_direction / n_direction_per_cone)
    for i in range(n_cone):
        coneProb = np.zeros(meshA.vertices.shape[0],dtype=float)
        for j in range(n_direction_per_cone):
            k = i*n_direction_per_cone+j
            vertex_function = np.dot(meshA.vertices,directions[k])
            if standardized:
                radius = np.linspace(-ball_radius,ball_radius,n_filtration)
            else:
                radius = np.linspace(np.amin(vertex_function),np.amax(vertex_function),n_filtration)        
            filtration = np.digitize(vertex_function,radius)
            outProb += pip[k*n_filtration+filtration]
    outProb /= n_direction
    return outProb      

def reconstruct_by_threshold(meshA, directions, pip, max_threshold = 0.1, n_cuts = 10, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0):
    n_direction = directions.shape[0]
    n_cone = int(n_direction/n_direction_per_cone)
    n_vertex = meshA.n_vertex
    vert_matrix = np.zeros(n_vertex,dtype=float)
    threshold_list = np.linspace(0,max_threshold,n_cuts)
    birth_time = np.linspace(0,1,n_cuts)
    for i_thr, threshold in enumerate(threshold_list):
        all_picked_vertices = []
        for i in range(n_cone):
            picked_vertices = np.zeros((n_direction_per_cone,n_vertex),dtype=bool)
            for j in range(n_direction_per_cone):
                k = i*n_direction_per_cone+j
                picked_filtrations = np.argwhere(pip[k*n_filtration:(k+1)*n_filtration] > threshold)
                vertex_function = np.dot(meshA.vertices,directions[k])
                radius = np.linspace(-ball_radius,ball_radius,n_filtration)
                filtration = np.digitize(vertex_function,radius)
                picked_vertices[j] = np.isin(filtration,picked_filtrations)
            all_picked_vertices = np.append(all_picked_vertices,np.argwhere(np.sum(picked_vertices,axis=0) >= 2))
        final_picked_vertices = np.unique(all_picked_vertices).astype(int)
        vert_matrix[final_picked_vertices] = birth_time[i_thr]
    return vert_matrix

def reconstruct_by_sorted_threshold(meshA, directions, pip, n_cuts = 10, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0, by_rank = False):
    n_direction = directions.shape[0]
    n_cone = int(n_direction/n_direction_per_cone)
    n_vertex = meshA.n_vertex
    vert_matrix = np.zeros(n_vertex,dtype=float)
    sorted_pip = np.sort(pip)
    n_pip = pip.size
    threshold_list = sorted_pip[::int(n_pip/n_cuts)]
    if by_rank:
        birth_time = np.linspace(0,1,n_cuts)
    else:
        birth_time = threshold_list / threshold_list[-1]
    for i_thr, threshold in enumerate(threshold_list):
        all_picked_vertices = []
        for i in range(n_cone):
            picked_vertices = np.zeros((n_direction_per_cone,n_vertex),dtype=bool)
            for j in range(n_direction_per_cone):
                k = i*n_direction_per_cone+j
                picked_filtrations = np.argwhere(pip[k*n_filtration:(k+1)*n_filtration] > threshold)
                vertex_function = np.dot(meshA.vertices,directions[k])
                radius = np.linspace(-ball_radius,ball_radius,n_filtration)
                filtration = np.digitize(vertex_function,radius)
                picked_vertices[j] = np.isin(filtration,picked_filtrations)
            all_picked_vertices = np.append(all_picked_vertices,np.argwhere(np.sum(picked_vertices,axis=0) >= n_direction_per_cone))
        final_picked_vertices = np.unique(all_picked_vertices).astype(int)
        if by_rank:
            vert_matrix[final_picked_vertices] = birth_time[i_thr]
        else:
            vert_matrix[final_picked_vertices] = threshold
    return vert_matrix

from scipy.stats import rankdata
def reconstruct_by_sorted_threshold_new(meshA, directions, pip, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0, by_rank = False):
    n_direction = directions.shape[0]
    n_cone = int(n_direction/n_direction_per_cone)
    n_vertex = meshA.n_vertex
    n_pip = pip.size
    pip_vert = np.zeros((n_vertex,n_cone,n_direction_per_cone),dtype=float)
    for i in range(n_cone):
        for j in range(n_direction_per_cone):
            k = i*n_direction_per_cone+j
            vertex_function = np.dot(meshA.vertices,directions[k])
            radius = np.linspace(-ball_radius,ball_radius,n_filtration)
            filtration = np.digitize(vertex_function,radius)
            pip_vert[:,i,j] = pip[k*n_filtration+filtration]
    height = np.amax(np.amin(pip_vert[:,:,:],axis=2),axis=1)
    if by_rank:
        rank = rankdata(height,method='dense')
        rank = rank/np.amax(rank)
        return rank
    else:
        return height

def reconstruct_from_lasso(meshA, directions, lasso_weights, n_cuts = 10, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0):
    n_direction = directions.shape[0]
    n_cone = int(n_direction/n_direction_per_cone)
    n_vertex = meshA.n_vertex
    vert_matrix = np.zeros(n_vertex,dtype=float)
    for threshold in np.linspace(0,1,n_cuts+1):
        all_picked_vertices = []
        for i in range(n_cone):
            picked_vertices = np.zeros((n_direction_per_cone,n_vertex),dtype=bool)
            for j in range(n_direction_per_cone):
                k = i*n_direction_per_cone+j
                lasso_indices = np.argwhere(lasso_weights > threshold)
                vertex_function = np.dot(meshA.vertices,directions[k])
                radius = np.linspace(-ball_radius,ball_radius,n_filtration)
                filtration = np.digitize(vertex_function,radius) + k*n_filtration
                picked_vertices[j] = np.isin(filtration,lasso_indices)
            all_picked_vertices = np.append(all_picked_vertices,np.argwhere(np.sum(picked_vertices,axis=0) > 1))
        final_picked_vertices = np.unique(all_picked_vertices).astype(int)
        vert_matrix[final_picked_vertices] = threshold
    return vert_matrix


