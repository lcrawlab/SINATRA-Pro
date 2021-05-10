#!/bin/python3
from traj_reader import *
from euler import *
from gp import *
from reconstruction import *
import sys

##########################################################################

# Name the variants, just for filename purpose
protA = "WT"
protB = "R164S"

# number of structures to draw from trajectory
n_sample = 1000

## Input files
struct_file_A = 'data/%s/md_0_1.gro'%protA
traj_file_A = 'data/%s/md_0_1_noPBC.xtc'%protA

struct_file_B = 'data/%s/md_0_1.gro'%protB
traj_file_B = 'data/%s/md_0_1_noPBC.xtc'%protB

## Simplicies construction parameters
## try using 2.0, 4.0 or 6.0 Angstrom, longer r takes longer computational time, so try smaller first, then larger to see if there are any differences
selection = 'protein'# and not (resid 164 and not backbone)'
sm_radius = 2.0      

## Filtration directions parameters
n_cone = 20
n_direction_per_cone = 8
cap_radius = 0.80

## Euler Characteristics (EC) calculation parameters
ec_type = "DECT"
n_filtration = 120

## Variable selection parameters
bandwidth = 0.01
#sampling_method = "ESS"
directory = "WT_R164S_whole"

## Parallelization setting
parallel = True  ## using multiple CPU cores for calculation
n_core = -1      ## Number of cores used, -1 = automatically detect and use all cores

## verbose
verbose = True

##########################################################################

## Read trajectory file and output aligned protein structures in pdb format
convert_traj_pdb_aligned(protA, protB, 
        struct_file_A=struct_file_A, 
        traj_file_A=traj_file_A, 
        struct_file_B=struct_file_B, 
        traj_file_B=traj_file_B, 
        align_frame=0, 
        n_sample=n_sample, 
        selection=selection, 
        offset=0, 
        directory=directory, 
        single=True) ## single="True" is for single run purpose, "False" for duplicate runs purpose which groups and names file with the frame offset.

#####################
### IF you already have your own aligned structure, start HERE put them in directory_pdb_A and directory_pdb_B, and put in the filename for the structure to use for visualization
#####################
directory_pdb_A = "%s/pdb/%s/"%(directory,protA)
directory_pdb_B = "%s/pdb/%s/"%(directory,protB)
reference_pdb_file = "%s/%s_frame0.pdb"%(directory_pdb_A,protA,protA) ## which pdb to use for visualization

## Converted protein structures into simplicial mesehes
convert_pdb_mesh(protA,protB,
        n_sample=n_sample, 
        sm_radius=sm_radius, 
        directory_pdb_A=directory_pdb_A, 
        directory_pdb_B=directory_pdb_B, 
        directory_mesh="%s/msh/"%(directory), 
        parallel=parallel, 
        n_core=n_core, 
        verbose=verbose)

## Calculate distributed cones of directions for EC calculations
directions = generate_equidistributed_cones(n_cone=n_cone,
        n_direction_per_cone=n_direction_per_cone,
        cap_radius=cap_radius,
        hemisphere=False)
#np.savetxt("%s/directions_%d_%d_%.2f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius),directions)

## EC calculations to convert simplicial meshes to topological summary statistics
X, y, not_vacuum = compute_ec_curve_folder(protA,protB,directions,
        n_sample=n_sample,
        ec_type=ec_type,
        n_cone=n_cone,
        n_direction_per_cone=n_direction_per_cone,
        cap_radius=cap_radius,
        n_filtration=n_filtration,
        sm_radius=sm_radius,
        directory_mesh_A = "%s/msh/%s_%.1f"%(directory,protA,sm_radius), 
        directory_mesh_B = "%s/msh/%s_%.1f"%(directory,protB,sm_radius), 
        parallel=parallel, 
        n_core=n_core, 
        verbose=verbose)

np.savetxt("%s/%s_%s_%s_%d_%d_%.2f_%d_norm_all.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration),X)
np.savetxt("%s/notvacuum_%s_%s_%s_%d_%d_%.2f_%d_norm_all.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration),not_vacuum)
np.savetxt('%s/%s_%s_label_all.txt'%(directory,protA,protB),y)    

## RATE calculation for variable selections from the topological summary statistics
kld, rates, delta, eff_samp_size = find_rate_variables_with_other_sampling_methods(X,y,
        bandwidth=bandwidth,
        parallel=parallel,
        n_core=n_core,
        verbose=verbose)
np.savetxt("%s/rates_%s_%s_%s_%d_%d_%.2f_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration),rates)

## reconstruct the RATE values onto the protein structures for visualization
## reconstruct probabilities are stored in "Temperature factor" column in the pdb format
## can then be visualized using Chimera or Pymol  
vert_prob = reconstruct_on_multiple_mesh(protA,protB,directions,
        rates=rates,
        not_vacuum=not_vacuum,
        n_sample=n_sample,
        n_direction_per_cone=n_direction_per_cone,
        n_filtration=n_filtration,
        sm_radius=sm_radius,
        directory_mesh="%s/msh/%s_%.1f"%(directory,protA,sm_radius),verbose=True)
np.savetxt("%s/vert_prob_DECT_%s_%s_%d_%d_%.2f_%d.txt"%(directory,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration),vert_prob)

write_vert_prob_on_pdb(vert_prob,
        protA=protA,
        protB=protB,
        selection=selection, 
        pdb_in_file=reference_pdb_file, 
        pdb_out_file="%s/vert_prob_DECT_%s_%s_%d_%d_%.2f_%d_all.pdb"%(directory,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration))
