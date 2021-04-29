#!/bin/python3

from mesh import *
import os, sys
import numpy as np
from scipy.special import logsumexp

# Reconstruction algorithms

from scipy.stats import rankdata
def reconstruct_by_sorted_threshold_new(meshfile, directions, rates, n_filtration = 25, n_direction_per_cone = 1, ball_radius = 1.0, by_rank = False, verbose = False):
    if verbose:
        sys.stdout.write('Reconstructing for %s ...\r'%meshfile)
        sys.stdout.flush()

    meshA = mesh()
    meshA.read_mesh_file(filename=meshfile)

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

def project_rate_on_nonvacuum(rates,not_vacuum):
    rates_new = np.zeros(not_vacuum.size,dtype=float)
    j = 0
    for i in range(not_vacuum.size):
        if not_vacuum[i]:
            rates_new[i] = rates[j]
            j += 1
    return rates_new

def reconstruct_on_multiple_mesh(protA, protB, directions, rates, not_vacuum, n_sample = 101, n_direction_per_cone = 4, n_filtration = 25, ball_radius = 1.0, directory_mesh = None, sm_radius = 4.0, by_rank = False, parallel = False, n_core = -1, verbose = False):
    if parallel:
        import multiprocessing
        from joblib import Parallel, delayed
        if n_core == -1:    
            n_core = multiprocessing.cpu_count()

    rates = project_rate_on_nonvacuum(rates,not_vacuum)
    if directory_mesh == None:
        directory_mesh = "%s_%s/mesh"%(protA,protB)
        if parallel:
            processed_list = Parallel(n_jobs=n_core)(delayed(reconstruct_by_sorted_threshold_new)('%s/%s_%.1f/%s_frame%d.msh'%(directory_mesh,protA,sm_radius,protA,frame), directions, rates, n_filtration, n_direction_per_cone, ball_radius, by_rank, verbose) for frame in range(n_sample))
            out_prob = np.array(processed_list)
        else:
            meshProtein = mesh()
            meshProtein.read_mesh_file(filename='%s/%s_%.1f/%s_frame0.msh'%(directory_mesh,protA,sm_radius,protA))
            out_prob = np.zeros((n_sample,meshProtein.vertices.shape[0]),dtype=float)
            for frame in range(n_sample):
                filename='%s/%s_%.1f/%s_frame%d.msh'%(directory_mesh,protA,sm_radius,protA,frame)
                out_prob[frame,:] = reconstruct_by_sorted_threshold_new('%s/%s_%.1f/%s_frame%d.msh'%(directory_mesh,protA,sm_radius,protA,frame), directions, rates, n_filtration = n_filtration, n_direction_per_cone = n_direction_per_cone, ball_radius = ball_radius, by_rank = by_rank, verbose = verbose)
        average_prob = np.average(out_prob,axis=0)
    else:
        if parallel:
            processed_list = Parallel(n_jobs=n_core)(delayed(reconstruct_by_sorted_threshold_new)(directory_mesh + '/' + filename, directions, rates, n_filtration, n_direction_per_cone, ball_radius, by_rank, verbose) for filename in os.listdir(directory_mesh))
            out_prob = np.array(processed_list)
        else:
            out_prob = []
            for filename in os.listdir(directory_mesh):
                if filename.endswith(".msh"):
                    prob = reconstruct_by_sorted_threshold_new(directory_mesh + '/' + filename, directions, rates, n_filtration = n_filtration, n_direction_per_cone = n_direction_per_cone, ball_radius = ball_radius, by_rank = by_rank, verbose = verbose)
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


