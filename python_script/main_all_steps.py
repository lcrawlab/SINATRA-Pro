#!/bin/python3

from traj_reader import *
from euler import *
from gp import *
from reconstruction import *

protA = "WT"
protB = "R164S"
nsample = 101

struct_file_A = '%s/md_0_1.gro'%protA
traj_file_A = '%s/md_0_1_noPBC.xtc'%protA

struct_file_B = '%s/md_0_1.gro'%protB
traj_file_B = '%s/md_0_1_noPBC.xtc'%protB

selection = 'resid 65:213'
sm_radius = 4.0
ec_type = "DECT"
n_filtration = 25
n_sample = 101

convert_traj_pdb_aligned(protA, protB, struct_file_A = struct_file_A, traj_file_A = traj_file_A, struct_file_B = struct_file_B, traj_file_B = traj_file_B, nsample = 101,selection=selection)

convert_pdb_mesh(protA,protB, n_sample = n_sample, radius = sm_radius)

directions = generate_equidistributed_cones(n_cone = 15, directions_per_cone = 4, cap_radius = 0.1, hemisphere = False)

X, y, notvacuum = compute_ec_curve_folder(protA,protB,directions,nsample=nsample,ec_type="ECT",n_cone=15,directions_per_cone=4,cap_radius=0.1,n_filtration=n_filtration,ball_radius=1.0,include_faces=True,sm_radius=radius)

kld, rates, delta, eff_samp_size = find_rate_variables_with_other_sampling_methods(X,y,bandwidth=0.01,sampling_method="ESS")

vert_prob = reconstruct_on_multiple_meshes(protA,protB,directions,rates,not_vacuum,n_sample=n_sample,directions=4,n_filtration=25,ball_radius=1.0,sm_radius=sm_radius)

write_vert_prob_on_pdb(vert_prob, pdb_in_file, pdb_out_file, selection = "protein")

