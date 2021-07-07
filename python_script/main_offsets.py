#!/bin/python3

from traj_reader import *
from euler import *
from gp import *
from reconstruction import *
import sys

protA = "WT"
protB = "R164S"
n_sample = 100

## Input files
struct_file_A = 'data/%s/md_0_1.gro'%protA
traj_file_A = 'data/%s/md_0_1_noPBC.xtc'%protA

struct_file_B = 'data/%s/md_0_1.gro'%protB
traj_file_B = 'data/%s/md_0_1_noPBC.xtc'%protB

## Simplicies construction parameters
#selection = 'protein'# and not (resid 164 and not backbone)'
selection = 'protein and resid 65:230'
sm_radius = 2.0

## Filtration directions parameters
n_cone = 20
n_direction_per_cone = 8
cap_radius = 0.80

## Euler Characteristics (EC) calculation parameters
ec_type = "DECT"
n_filtration = 120
parallel = True ## using multiple CPU cores for calculation
n_core = -1 ## Number of cores used for EC calculation

## Variable selection parameters
bandwidth = 0.01
sampling_method = "ESS"
directory = "WT_R164S_65_230"

n_step = int(10000/n_sample)
n_skip = int(n_step/10)

for offset in range(0,n_step,n_skip):
        
    #print(offset)    
    print("%d %d %.2f %d %d"%(n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))   
    
    if os.path.exists("%s/vert_prob_DECT_%s_%s_%d_%d_%.2f_%d_offset_%d.txt"%(directory,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset)):
        continue

    ## Read trajectory file and output aligned protein structures in pdb format
    convert_traj_pdb_aligned(protA, protB, struct_file_A = struct_file_A, traj_file_A = traj_file_A, struct_file_B = struct_file_B, traj_file_B = traj_file_B, align_frame = offset, n_sample = n_sample,selection = selection, offset = offset, directory = directory)
    #convert_traj_pdb_aligned(protA, protB, struct_file_A = struct_file_A, traj_file_A = traj_file_A, struct_file_B = struct_file_B, traj_file_B = traj_file_B, align_frame = 10000-100+offset, n_sample = n_sample,selection = selection, offset = offset, directory = directory)

    ## Converted protein structures into simplicial mesehes
    convert_pdb_mesh(protA, protB, n_sample = n_sample, sm_radius = sm_radius, directory_pdb_A = "%s/pdb/%s_offset_%d"%(directory,protA,offset), directory_pdb_B = "%s/pdb/%s_offset_%d"%(directory,protB,offset), directory_mesh = "%s/msh_offset_%d/"%(directory,offset), parallel = True, n_core = n_core)
    
    ## Calculate distributed cones of directions for EC calculations
    directions = generate_equidistributed_cones(n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,hemisphere=True)
    #np.savetxt("%s/directions_%d_%d_%.2f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius),directions)

    if os.path.exists("%s/%s_%s_%s_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset)):
        X = np.loadtxt("%s/%s_%s_%s_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
        not_vacuum = np.loadtxt("%s/notvacuum_%s_%s_%s_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
        y = np.loadtxt('%s/%s_%s_label.txt'%(directory,protA,protB))
    else:
        ## EC calculations to convert simplicial meshes to topological summary statistics
        X, y, not_vacuum = compute_ec_curve_folder(protA,protB,directions,n_sample=n_sample,ec_type=ec_type,n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,n_filtration=n_filtration,sm_radius=sm_radius,directory_mesh_A = "%s/msh_offset_%d/%s_%.1f"%(directory,offset,protA,sm_radius), directory_mesh_B = "%s/msh_offset_%d/%s_%.1f"%(directory,offset,protB,sm_radius), parallel=parallel,n_core=n_core)
        np.savetxt("%s/%s_%s_%s_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),X)
        np.savetxt("%s/notvacuum_%s_%s_%s_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),not_vacuum)
        np.savetxt('%s/%s_%s_label.txt'%(directory,protA,protB),y)    
    
    if os.path.exists("%s/rates_%s_%s_%s_%d_%d_%.2f_%d_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset)):
        rates = np.loadtxt("%s/rates_%s_%s_%s_%d_%d_%.2f_%d_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
    else:
        ## RATE calculation for variable selections from the topological summary statistics
        kld, rates, delta, eff_samp_size = find_rate_variables_with_other_sampling_methods(X,y,bandwidth=bandwidth,sampling_method=sampling_method,parallel=True,n_core=n_core)
        np.savetxt("%s/rates_%s_%s_%s_%d_%d_%.2f_%d_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),rates)
    
    if not os.path.exists("%s/vert_prob_DECT_%s_%s_%d_%d_%.2f_%d_offset_%d.txt"%(directory,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset)):
        ## reconstruct the RATE values onto the protein structures for visualization
        ## reconstruct probabilities are stored in "Temperature factor" column in the pdb format
        ## can then be visualized using Chimera or Pymol  
        vert_prob = reconstruct_on_multiple_mesh(protA,protB,directions,rates,not_vacuum,n_sample=n_sample,n_direction_per_cone=n_direction_per_cone,n_filtration=n_filtration,sm_radius=sm_radius,directory_mesh="%s/msh_offset_%d/%s_%.1f"%(directory,offset,protA,sm_radius))
        np.savetxt("%s/vert_prob_DECT_%s_%s_%d_%d_%.2f_%d_offset_%d.txt"%(directory,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),vert_prob)
        write_vert_prob_on_pdb(vert_prob,protA=protA,protB=protB,selection=selection,pdb_in_file = "%s/pdb/%s_offset_%d/%s_frame0.pdb"%(directory,protA,offset,protA),pdb_out_file="%s/vert_prob_DECT_%s_%s_%d_%d_%.2f_%d_offset_%d.pdb"%(directory,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
