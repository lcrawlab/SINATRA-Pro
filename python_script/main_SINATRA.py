#!/bin/python3

from traj_reader import *
from euler import *
from gp import *
from reconstruction import *

protA = "WT"
protB = "R164S"
nsample = 101

## Input files
struct_file_A = 'data/%s/md_0_1.gro'%protA
traj_file_A = 'data/%s/md_0_1_noPBC.xtc'%protA

struct_file_B = 'data/%s/md_0_1.gro'%protB
traj_file_B = 'data/%s/md_0_1_noPBC.xtc'%protB

## Simplicies construction parameters
selection = 'resid 65:230'
sm_radius = 4.0

## Filtration directions parameters
n_cone = 15
n_direction_per_cone = 8
cap_radius = 0.1

## Euler Characteristics (EC) calculation parameters
ec_type = "DECT"
n_filtration = 55
n_sample = 100
parallel = True ## using multiple CPU cores for calculation
n_core = -1 ## Number of cores used for EC calculation

## Variable selection parameters
bandwidth = 0.01
sampling_method = "ESS"

for offset in range(0,100,10):

    ## Read trajectory file and output aligned protein structures in pdb format
    convert_traj_pdb_aligned(protA, protB, struct_file_A = struct_file_A, traj_file_A = traj_file_A, struct_file_B = struct_file_B, traj_file_B = traj_file_B, align_frame = offset, n_sample = n_sample,selection = selection, offset = offset)

    ## Converted protein structures into simplicial mesehes
    convert_pdb_mesh(protA, protB, n_sample = n_sample, sm_radius = sm_radius, directory_pdb_A = "%s_%s/pdb/%s_offset_%d"%(protA,protB,protA,offset), directory_pdb_B = "%s_%s/pdb/%s_offset_%d"%(protA,protB,protB,offset), directory_mesh = "%s_%s/msh_offset_%d/"%(protA,protB,offset), parallel = True, n_core = n_core)

    ## Calculate distributed cones of directions for EC calculations
    directions = generate_equidistributed_cones(n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,hemisphere=True)
    np.savetxt("WT_R164S/directions_%d_%d_0.1.txt"%(n_cone,n_direction_per_cone),directions)

    ## EC calculations to convert simplicial meshes to topological summary statistics
    X, y, not_vacuum = compute_ec_curve_folder(protA,protB,directions,n_sample=n_sample,ec_type=ec_type,n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,n_filtration=n_filtration,sm_radius=sm_radius,directory_mesh_A = "%s_%s/msh_offset_%d/%s_%.1f"%(protA,protB,offset,protA,sm_radius), directory_mesh_B = "%s_%s/msh_offset_%d/%s_%.1f"%(protA,protB,offset,protB,sm_radius), parallel=parallel,n_core=n_core)

    np.savetxt("%s_%s/%s_%s_%s_%d_%d_%.1f_%d_norm_offset_%d.txt"%(protA,protB,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),X)
    np.savetxt("%s_%s/notvacuum_%s_%s_%s_%d_%d_%.1f_%d_norm_offset_%d.txt"%(protA,protB,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),not_vacuum)
    np.savetxt('%s_%s/%s_%s_label.txt'%(protA,protB,protA,protB),y)    

    #X = np.loadtxt("WT_R164S/DECT_WT_R164S_15_4_0.1_25_norm.txt")
    #not_vacuum = np.loadtxt("WT_R164S/notvacuum_DECT_WT_R164S_15_4_0.1_25_norm.txt")
    #y = np.loadtxt('WT_R164S/WT_R164S_label.txt')

    ## RATE calculation for variable selections from the topological summary statistics
    kld, rates, delta, eff_samp_size = find_rate_variables_with_other_sampling_methods(X,y,bandwidth=bandwidth,sampling_method=sampling_method,parallel=True,n_core=n_core)
    np.savetxt("%s_%s/rates_%s_%s_%s_%d_%d_0.1_%d_offset_%d.txt"%(protA,protB,ec_type,protA,protB,n_cone,n_direction_per_cone,n_filtration,offset),rates)

    ## reconstruct the RATE values onto the protein structures for visualization
    ## reconstruct probabilities are stored in "Temperature factor" column in the pdb format
    ## can then be visualized using Chimera or pymol 
    rates = np.loadtxt("%s_%s/rates_DECT_%s_%s_%d_%d_0.1_%d_offset_%d.txt"%(protA,protB,protA,protB,n_cone,n_direction_per_cone,n_filtration,offset))
    directions = np.loadtxt("WT_R164S/directions_%d_%d_0.1.txt"%(n_cone,n_direction_per_cone))
    vert_prob = reconstruct_on_multiple_mesh(protA,protB,directions,rates,not_vacuum,n_sample=n_sample,n_direction_per_cone=n_direction_per_cone,n_filtration=n_filtration,sm_radius=sm_radius,directory_mesh="%s_%s/msh_offset_%d/%s_%.1f"%(protA,protB,offset,protA,sm_radius))
    np.savetxt("%s_%s/vert_prob_DECT_%s_%s_%d_%d_%.1f_%d_offset_%d.txt"%(protA,protB,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),vert_prob)
    write_vert_prob_on_pdb(vert_prob,protA=protA,protB=protB,selection=selection,pdb_in_file = "%s_%s/pdb/%s_offset_%d/%s_frame0.pdb"%(protA,protB,protA,offset,protA),pdb_out_file="%s_%s/vert_prob_DECT_%s_%s_%d_%d_%.1f_%d_offset_%d.pdb"%(protA,protB,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
