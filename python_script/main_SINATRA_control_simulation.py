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
n_direction_per_cone = int(sys.argv[1])
cap_radius = 0.1

## Euler Characteristics (EC) calculation parameters
ec_type = "DECT"
n_filtration = 25
n_sample = 100
parallel = True ## using multiple CPU cores for calculation
n_core = -1 ## Number of cores used for EC calculation

## Variable selection parameters
bandwidth = 0.01
sampling_method = "ESS"

import multiprocessing
n_core = multiprocessing.cpu_count()
print("Detected %d cores"%n_core)

## Input files
struct_file = 'data/WT/md_0_1.gro'
traj_file = 'data/WT/md_0_1_noPBC.xtc'

directory = "simulation_removed_omega_loop_offset"
if not os.path.exists(directory):
    os.mkdir(directory)

for offset in range(0,100,10):

    ## Pick different set of frames from trajectory
    #convert_traj_pdb_single_offset(prot, struct_file, traj_file, align_frame = offset, n_sample = 101, offset = offset, selection = selection, directory = directory)

    ## Perturb pdb positions and convert to mesh
    #remove_protein_region(sm_radius = 4.0, selection='resid 163:178', prot = prot, n_sample = n_sample, directory_original = "%s/pdb/offset_%d"%(directory,offset), directory_pdb_B = "%s/pdb/offset_%d_perturbed"%(directory,offset), directory_new = directory, directory_mesh_A = "%s/msh/offset_%d_original"%(directory,offset), directory_mesh_B = "%s/msh/offset_%d_perturbed"%(directory,offset), parallel = True, n_core = n_core)

    ## Calculate distributed cones of directions for EC calculations
    directions = generate_equidistributed_cones(n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,hemisphere=False)
    np.savetxt("%s/directions_%d_%d_%.1f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius),directions)

    #directions = np.loadtxt("%s/directions_%d_%d_%.1f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius))

    ## EC calculations to convert simplicial meshes to topological summary statistics
    X, y, not_vacuum = compute_ec_curve_folder(protA,protB,directions,directory_mesh_A="%s/msh/offset_%d_original"%(directory,offset), directory_mesh_B="%s/msh/offset_%d_perturbed"%(directory,offset),n_sample=n_sample,ec_type=ec_type,n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,n_filtration=n_filtration,sm_radius=sm_radius,parallel=parallel,n_core=n_core)

    np.savetxt("%s/DECT_original_perturbed_%d_%d_%.1f_%d_norm_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),X)
    np.savetxt('%s/original_perturbed_label.txt'%directory,y)
    np.savetxt("%s/notvacuum_DECT_original_perturbed_%d_%d_%.1f_%d_norm_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),not_vacuum)

    #X = np.loadtxt("%s/DECT_original_perturbed_%d_%d_%.1f_%d_norm.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration))
    #y = np.loadtxt('%s/original_perturbed_label.txt'%directory)
    #not_vacuum = np.loadtxt("%s/notvacuum_DECT_original_perturbed_%d_%d_%.1f_%d_norm.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration))

    ## RATE calculation for variable selections from the topological summary statistics
    kld, rates, delta, eff_samp_size = find_rate_variables_with_other_sampling_methods(X,y,bandwidth=bandwidth,sampling_method=sampling_method,parallel=True,n_core=n_core)
    np.savetxt("%s/rates_DECT_original_perturbed_%d_%d_%.1f_%d_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),rates)

    ## reconstruct the RATE values onto the protein structures for visualization
    ## reconstruct probabilities are stored in "Temperature factor" column in the pdb format
    ## can then be visualized using Chimera or pymol 
    rates = np.loadtxt("%s/rates_DECT_original_perturbed_%d_%d_%.1f_%d_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
    vert_prob = reconstruct_on_multiple_mesh(protA,protB,directions,rates,not_vacuum,n_sample=n_sample,directory_mesh="%s/msh/offset_%d_original"%(directory,offset),n_direction_per_cone=n_direction_per_cone,n_filtration=n_filtration,sm_radius=sm_radius)
    np.savetxt("%s/vert_prob_DECT_original_perturbed_%d_%d_%.1f_%d_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),vert_prob)

    vert_prob = np.loadtxt("%s/vert_prob_DECT_original_perturbed_%d_%d_%.1f_%d_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
    write_vert_prob_on_pdb(vert_prob,protA=protA,protB=protB,selection=selection,pdb_in_file="%s/pdb/offset_%d/%s_frame0.pdb"%(directory,offset,prot),pdb_out_file="%s/vert_prob_DECT_original_perturbed_%d_%d_%.1f_%d_offset_%d.pdb"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))

