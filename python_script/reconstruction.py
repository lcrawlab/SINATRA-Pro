#!/bin/python3

from mesh import *
import os, sys
import numpy as np
from scipy.special import logsumexp

# Reconstruction algorithms

def reconstruct_by_average(meshA, directions, rates, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0, standardized = True):
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
            outProb += rates[k*n_filtration+filtration]
    outProb /= n_direction
    return outProb      

def reconstruct_by_threshold(meshA, directions, rates, max_threshold = 0.1, n_cuts = 10, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0):
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
                picked_filtrations = np.argwhere(rates[k*n_filtration:(k+1)*n_filtration] > threshold)
                vertex_function = np.dot(meshA.vertices,directions[k])
                radius = np.linspace(-ball_radius,ball_radius,n_filtration)
                filtration = np.digitize(vertex_function,radius)
                picked_vertices[j] = np.isin(filtration,picked_filtrations)
            all_picked_vertices = np.append(all_picked_vertices,np.argwhere(np.sum(picked_vertices,axis=0) >= 2))
        final_picked_vertices = np.unique(all_picked_vertices).astype(int)
        vert_matrix[final_picked_vertices] = birth_time[i_thr]
    return vert_matrix

def reconstruct_by_sorted_threshold(meshA, directions, rates, n_cuts = 10, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0, by_rank = False):
    n_direction = directions.shape[0]
    n_cone = int(n_direction/n_direction_per_cone)
    n_vertex = meshA.n_vertex
    vert_matrix = np.zeros(n_vertex,dtype=float)
    sorted_rates = np.sort(rates)
    n_rates = rates.size
    threshold_list = sorted_rates[::int(n_rates/n_cuts)]
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
                picked_filtrations = np.argwhere(rates[k*n_filtration:(k+1)*n_filtration] > threshold)
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
def reconstruct_by_sorted_threshold_new(meshA, directions, rates, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0, by_rank = False):
    n_direction = directions.shape[0]
    n_cone = int(n_direction/n_direction_per_cone)
    n_vertex = meshA.n_vertex
    n_rates = rates.size
    rates_vert = np.zeros((n_vertex,n_cone,n_direction_per_cone),dtype=float)
    for i in range(n_cone):
        for j in range(n_direction_per_cone):
            k = i*n_direction_per_cone+j
            vertex_function = np.dot(meshA.vertices,directions[k])
            radius = np.linspace(-ball_radius,ball_radius,n_filtration)
            filtration = np.digitize(vertex_function,radius)-1
            rates_vert[:,i,j] = rates[k*n_filtration+filtration]
    height = np.amax(np.amin(rates_vert[:,:,:],axis=2),axis=1)
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


def project_rate_on_nonvacuum(rates,not_vacuum):
    rates_new = np.zeros(not_vacuum.size,dtype=float)
    j = 0
    for i in range(not_vacuum.size):
        if not_vacuum[i]:
            rates_new[i] = rates[j]
            j += 1
    return rates_new


def reconstruct_on_multiple_mesh(protA, protB, directions, rates, not_vacuum, n_sample = 101, n_direction_per_cone = 4, n_filtration = 25, ball_radius = 1.0, directory_mesh = None, sm_radius = 4.0, parallel = False, n_core = -1, verbose = False):
    rates = project_rate_on_nonvacuum(rates,not_vacuum)
    if directory_mesh == None:
        directory_mesh = "%s_%s/mesh"%(protA,protB)
        meshProtein = mesh()
        meshProtein.read_mesh_file(filename='%s/%s_%.1f/%s_frame0.msh'%(directory_mesh,protA,sm_radius,protA))
        out_prob = np.zeros((n_sample,meshProtein.vertices.shape[0]),dtype=float)
        for frame in range(n_sample):
            if verbose:
                sys.stdout.write('Reconstructing for Frame %d...\r'%frame)
                sys.stdout.flush()
            meshProtein = mesh()
            meshProtein.read_mesh_file(filename='%s/%s_%.1f/%s_frame%d.msh'%(directory_mesh,protA,sm_radius,protA,frame))
            out_prob[frame,:] = reconstruct_by_sorted_threshold_new(meshProtein, directions, rates, n_filtration = n_filtration, n_direction_per_cone = n_direction_per_cone, ball_radius = ball_radius, by_rank = False)
        average_prob = np.average(out_prob,axis=0)
    else:
        out_prob = []
        for filename in os.listdir(directory_mesh):
            if filename.endswith(".msh"):
                if verbose:
                    sys.stdout.write('Reconstructing for mesh %s...\r'%filename)
                    sys.stdout.flush()
                meshProtein = mesh()
                meshProtein.read_mesh_file(filename=directory_mesh + '/' + filename)
                prob = reconstruct_by_sorted_threshold_new(meshProtein, directions, rates, n_filtration = n_filtration, n_direction_per_cone = n_direction_per_cone, ball_radius = ball_radius, by_rank = False)
                out_prob.append(prob)
        out_prob = np.array(out_prob)
        average_prob = np.average(out_prob,axis=0)
    return average_prob

def write_vert_prob_on_pdb(vert_prob, protA = None, protB = None, pdb_in_file = None, pdb_out_file = None, selection = "protein"):
    import MDAnalysis as mda
    if selection == None:
        selection = "protein"
    if pdb_in_file == None:
        pdb_in_file = "%s_%s/pdb/%s/%s_frame0.pdb"%(protA,protB,protA,protA)
    if pdb_out_file == None:
        pdb_out_file = "%s_%s/%s_reconstructed.pdb"%(protA,protB,protA)
    u = mda.Universe(pdb_in_file)
    protein = u.select_atoms(selection)
    u.add_TopologyAttr('tempfactors')
    ymin = np.amin(vert_prob)
    ymax = np.amax(vert_prob)
    y = (vert_prob - ymin)/(ymax-ymin)*100
    protein.tempfactors = y
    protein.write(pdb_out_file)
    return


