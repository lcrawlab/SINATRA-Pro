#!/bin/python3

from traj_reader import *
from euler import *
from gp import *
from reconstruction import *
import sys

protA = "WT"
protB = "R164S"
n_sample = 1000

## Input files
struct_file_A = 'data/%s/md_0_1.gro'%protA
traj_file_A = 'data/%s/md_0_1_noPBC.xtc'%protA

struct_file_B = 'data/%s/md_0_1.gro'%protB
traj_file_B = 'data/%s/md_0_1_noPBC.xtc'%protB

## Simplicies construction parameters
selection = 'protein'# and not (resid 164 and not backbone)'
sm_radius = 2.0 #float(sys.argv[1])

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
directory = "%s_%s_whole_2.0"%(protA,protB)

## Read trajectory file and output aligned protein structures in pdb format
convert_traj_pdb_aligned(protA, protB, struct_file_A = struct_file_A, traj_file_A = traj_file_A, struct_file_B = struct_file_B, traj_file_B = traj_file_B, align_frame = 0, n_sample = n_sample,selection = selection, offset = 0, directory = directory, single = True)

## Converted protein structures into simplicial mesehes
convert_pdb_mesh(protA, protB, n_sample = n_sample, sm_radius = sm_radius, directory_pdb_A = "%s/pdb/%s/"%(directory,protA), directory_pdb_B = "%s/pdb/%s/"%(directory,protB), directory_mesh = "%s/msh_all/"%(directory), parallel = True, n_core = n_core)

## Calculate distributed cones of directions for EC calculations
directions = generate_equidistributed_cones(n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,hemisphere=True)
#np.savetxt("%s/directions_%d_%d_%.2f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius),directions)

## EC calculations to convert simplicial meshes to topological summary statistics
X, y, not_vacuum = compute_ec_curve_folder(protA,protB,directions,n_sample=n_sample,ec_type=ec_type,n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,n_filtration=n_filtration,sm_radius=sm_radius,directory_mesh_A = "%s/msh_all/%s_%.1f"%(directory,protA,sm_radius), directory_mesh_B = "%s/msh_all/%s_%.1f"%(directory,protB,sm_radius), parallel=parallel,n_core=n_core)
np.savetxt("%s/%s_%s_%s_%d_%d_%.2f_%d_norm_all.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration),X)
np.savetxt("%s/notvacuum_%s_%s_%s_%d_%d_%.2f_%d_norm_all.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration),not_vacuum)
np.savetxt('%s/%s_%s_label_all.txt'%(directory,protA,protB),y)    

## RATE calculation for variable selections from the topological summary statistics
kld, rates, delta, eff_samp_size = find_rate_variables_with_other_sampling_methods(X,y,bandwidth=bandwidth,sampling_method=sampling_method,parallel=True,n_core=n_core)
np.savetxt("%s/rates_%s_%s_%s_%d_%d_%.2f_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration),rates)

## reconstruct the RATE values onto the protein structures for visualization
## reconstruct probabilities are stored in "Temperature factor" column in the pdb format
## can then be visualized using Chimera or Pymol  
vert_prob = reconstruct_on_multiple_mesh(protA,protB,directions,rates,not_vacuum,n_sample=n_sample,n_direction_per_cone=n_direction_per_cone,n_filtration=n_filtration,sm_radius=sm_radius,directory_mesh="%s/msh_all/%s_%.1f"%(directory,protA,sm_radius))
np.savetxt("%s/vert_prob_DECT_%s_%s_%d_%d_%.2f_%d.txt"%(directory,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration),vert_prob)
write_vert_prob_on_pdb(vert_prob,protA=protA,protB=protB,selection=selection,pdb_in_file = "%s/pdb/%s/%s_frame0.pdb"%(directory,protA,protA),pdb_out_file="%s/vert_prob_DECT_%s_%s_%d_%d_%.2f_%d_all.pdb"%(directory,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration))
