#!/bin/python3

import os, sys
from traj_reader import *
from euler import *
from gp import *
from reconstruction import *
from control_simulation import *

prot = "WT"
nsample = 101
protA = "original"
protB = "perturbed"

## Simplicies construction parameters
selection = 'resid 65:213'
sm_radius = 4.0

## Filtration directions parameters
n_cone = 15
n_direction_per_cone = 8
cap_radius = 0.1

## Euler Characteristics (EC) calculation parameters
ec_type = "DECT"
n_filtration = 25
n_sample = 101
parallel = True ## using multiple CPU cores for calculation
n_core = -1 ## Number of cores used for EC calculation

## Variable selection parameters
bandwidth = 0.01
sampling_method = "ESS"

import multiprocessing
n_core = multiprocessing.cpu_count()
print("Detected %d cores"%n_core)

for i in np.arange(65,213,20):
    start = i
    end = i + 20
    directory = "simulation_removed_%d_%d"%(start,end)
    if not os.path.exists(directory):
        os.mkdir(directory)
        
    ## Perturb pdb positions and convert to mesh
    remove_protein_region(sm_radius = 4.0, selection='resid %d:%d'%(start,end), prot = prot, n_sample = n_sample, directory_original = "WT_R164S/pdb/WT/", directory_new = directory, parallel = True, n_core = n_core)

    ## Calculate distributed cones of directions for EC calculations
    directions = generate_equidistributed_cones(n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,hemisphere=False)
    np.savetxt("%s/directions_%d_%d_%.1f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius),directions)

    #directions = np.loadtxt("%s/directions_%d_%d_%.1f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius))

    ## EC calculations to convert simplicial meshes to topological summary statistics
    X, y, not_vacuum = compute_ec_curve_folder(protA,protB,directions,directory_mesh_A="%s/msh/original"%directory, directory_mesh_B="%s/msh/perturbed"%directory,n_sample=n_sample,ec_type=ec_type,n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,n_filtration=n_filtration,sm_radius=sm_radius,parallel=parallel,n_core=n_core)

    np.savetxt("%s/DECT_original_perturbed_%d_%d_%.1f_%d_norm.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration),X)
    np.savetxt('%s/original_perturbed_label.txt'%directory,y)
    np.savetxt("%s/notvacuum_DECT_original_perturbed_%d_%d_%.1f_%d_norm.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration),not_vacuum)

    #X = np.loadtxt("%s/DECT_original_perturbed_%d_%d_%.1f_%d_norm.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration))
    #y = np.loadtxt('%s/original_perturbed_label.txt'%directory)
    #not_vacuum = np.loadtxt("%s/notvacuum_DECT_original_perturbed_%d_%d_%.1f_%d_norm.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration))

    ## RATE calculation for variable selections from the topological summary statistics
    kld, rates, delta, eff_samp_size = find_rate_variables_with_other_sampling_methods(X,y,bandwidth=bandwidth,sampling_method=sampling_method,parallel=True,n_core=n_core)
    np.savetxt("%s/rates_DECT_original_perturbed_%d_%d_%.1f_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration),rates)

    ## reconstruct the RATE values onto the protein structures for visualization
    ## reconstruct probabilities are stored in "Temperature factor" column in the pdb format
    ## can then be visualized using Chimera or pymol 
    rates = np.loadtxt("%s/rates_DECT_original_perturbed_%d_%d_%.1f_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration))
    vert_prob = reconstruct_on_multiple_mesh(protA,protB,directions,rates,not_vacuum,n_sample=n_sample,directory_mesh="%s/msh/original"%directory,n_direction_per_cone=n_direction_per_cone,n_filtration=n_filtration,sm_radius=sm_radius)
    np.savetxt("%s/vert_prob_DECT_original_perturbed_%d_%d_%.1f_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration),vert_prob)

    vert_prob = np.loadtxt("%s/vert_prob_DECT_original_perturbed_%d_%d_%.1f_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration))
    write_vert_prob_on_pdb(vert_prob,protA=protA,protB=protB,selection=selection,pdb_in_file="WT_R164S/pdb/WT/WT_frame0.pdb",pdb_out_file="%s/vert_prob_DECT_original_perturbed_%d_%d_%.1f_%d.pdb"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration))

