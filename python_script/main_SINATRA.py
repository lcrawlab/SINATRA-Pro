#!/bin/python3

from traj_reader import *
from euler import *
from gp import *
from reconstruction import *

protA = "WT"
protB = "R164S"
nsample = 101

## Input files
struct_file_A = '../%s/md_0_1.gro'%protA
traj_file_A = '../%s/md_0_1_noPBC.xtc'%protA

struct_file_B = '../%s/md_0_1.gro'%protB
traj_file_B = '../%s/md_0_1_noPBC.xtc'%protB

## Simplicies construction parameters
selection = 'resid 65:230'
sm_radius = 4.0

## Filtration directions parameters
n_cone = 15
n_direction_per_cone = 4
cap_radius = 0.1

## Euler Characteristics (EC) calculation parameters
ec_type = "DECT"
n_filtration = 25
n_sample = 101
parallel = False ## using multiple CPU cores for calculation
n_core = 4 ## Number of cores used for EC calculation

## Variable selection parameters
bandwidth = 0.01
sampling_method = "ESS"

## Read trajectory file and output aligned protein structures in pdb format
#convert_traj_pdb_aligned(protA, protB, struct_file_A = struct_file_A, traj_file_A = traj_file_A, struct_file_B = struct_file_B, traj_file_B = traj_file_B, align_frame = 0, nsample = 101,selection = selection)

## Converted protein structures into simplicial mesehes
convert_pdb_mesh(protA,protB, n_sample = n_sample, sm_radius = sm_radius, parallel = True, n_core = 4)

## Calculate distributed cones of directions for EC calculations
directions = generate_equidistributed_cones(n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,hemisphere=False)
np.savetxt("WT_R164S/directions_4_15_0.1.txt",directions)

directions = np.loadtxt("WT_R164S/directions_4_15_0.1.txt")
## EC calculations to convert simplicial meshes to topological summary statistics
X, y, not_vacuum = compute_ec_curve_folder(protA,protB,directions,n_sample=n_sample,ec_type=ec_type,n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,n_filtration=n_filtration,sm_radius=sm_radius,parallel=parallel,n_core=n_core)

#X = np.loadtxt("WT_R164S/DECT_WT_R164S_15_4_0.1_25_norm.txt")
#not_vacuum = np.loadtxt("WT_R164S/notvacuum_DECT_WT_R164S_15_4_0.1_25_norm.txt")
#y = np.loadtxt('WT_R164S/WT_R164S_label.txt')

## RATE calculation for variable selections from the topological summary statistics
kld, rates, delta, eff_samp_size = find_rate_variables_with_other_sampling_methods(X,y,bandwidth=bandwidth,sampling_method=sampling_method)
np.savetxt("WT_R164S/rates_DECT_WT_R164S_15_4_0.1_25.txt",rates)

## reconstruct the RATE values onto the protein structures for visualization
## reconstruct probabilities are stored in "Temperature factor" column in the pdb format
## can then be visualized using Chimera or pymol 
#rates = np.loadtxt("WT_R164S/rates_DECT_WT_R164S_15_4_0.1_25.txt")
#directions = np.loadtxt("WT_R164S/directions_4_15_0.1.txt")
vert_prob = reconstruct_on_multiple_mesh(protA,protB,directions,rates,not_vacuum,n_sample=n_sample,n_direction_per_cone=4,n_filtration=25,sm_radius=sm_radius)
write_vert_prob_on_pdb(vert_prob,protA=protA,protB=protB,selection=selection)

