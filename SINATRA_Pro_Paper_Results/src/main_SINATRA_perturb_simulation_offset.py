#!/bin/python3

import os, sys
from traj_reader import *
from euler import *
from gp import *
from reconstruction import *
from control_simulation import *

prot = "WT"
protA = "original"
protB = "perturbed"

## Simplicies construction parameters
selection = 'resid 65:213'
sm_radius = 4.0

#dr = 0.0
#dr_vector = [0.5,0.5,0.5]
#atom_noise = 0.0
#perturb_noise = 0.0

r = 0.5
noise = 0.0

## Filtration directions parameters
#n_cone = 20
#n_direction_per_cone = 8
#cap_radius = 0.30

## Euler Characteristics (EC) calculation parameters
ec_type = "DECT"
#n_filtration = 60
n_sample = 100
parallel = True ## using multiple CPU cores for calculation
n_core = -1 ## Number of cores used for EC calculation

#n_cone = 20
#n_direction_per_cone = 8
#cap_radius = 0.30
#n_filtration = 60
n_cone = 20
n_direction_per_cone = 8
cap_radius = 0.80
n_filtration = 120

## Variable selection parameters
bandwidth = 0.01
sampling_method = "ESS"

import multiprocessing
n_core = multiprocessing.cpu_count()
print("Detected %d cores"%n_core)

## Input files
struct_file = 'data/WT/md_0_1.gro'
traj_file = 'data/WT/md_0_1_noPBC.xtc'

#directory_original = "WT_R164S_65_213"
#directory_perturb = "simulation_perturb_omega_loop_vector_%.1f_%.1f_%.1f_%.1f_offset"%(atom_noise,dr_vector[0],dr_vector[1],dr_vector[2])
#directory = "simulation_perturb_omega_loop_vector_%.1f_%.1f_%.1f_%.1f_offset"%(atom_noise,dr_vector[0],dr_vector[1],dr_vector[2])

directory_perturb = "simulation_perturb_omega_loop_sphere_%.1f_%.1f_offset_0.4ns"%(r,noise)
directory = directory_perturb #"simulation_perturb_omega_loop_sphere_%.1f_%.1f_offset"%(r,noise)

if not os.path.exists(directory_perturb):
    os.mkdir(directory_perturb)
if not os.path.exists(directory_perturb + '/pdb'):
    os.mkdir(directory_perturb + '/pdb')
if not os.path.exists(directory_perturb + '/msh'):
    os.mkdir(directory_perturb + '/msh')

if not os.path.exists(directory):
    os.mkdir(directory)

if not os.path.exists(directory):
    os.mkdir(directory)

step = int(10000/n_sample)
skip = int(step/10)


for offset in range(0,50,10):
    ## Original set of frames
    convert_traj_pdb_aligned(protA = prot, protB = prot, struct_file_A = struct_file, traj_file_A = traj_file, struct_file_B = struct_file, traj_file_B = traj_file, align_frame = offset, n_sample = n_sample, offset = offset, selection = selection, directory = directory)
    ## Pick different set of frames from trajectory
    convert_traj_pdb_aligned(protA = prot, protB = prot, struct_file_A = struct_file, traj_file_A = traj_file, struct_file_B = struct_file, traj_file_B = traj_file, align_frame = offset, n_sample = n_sample, offset = offset+40, selection = selection, directory = directory)   
    ## Perturb pdb positions
    #perturb_protein_region_both(sm_radius = 4.0, selection='resid 163:178', atom_noise = atom_noise, dr = dr, dr_vector = dr_vector, prot = prot, n_sample = n_sample, directory_original = "%s/pdb/WT_offset_%d"%(directory,offset+50), directory_pdb_B = "%s/pdb/offset_%d_perturbed"%(directory,offset+50), directory_new = directory, parallel = True, n_core = n_core)
    perturb_protein_region_sphere(sm_radius = 4.0, selection='resid 163:178', noise = noise, r = r, prot = prot, n_sample = n_sample, directory_original = "%s/pdb/WT_offset_%d"%(directory,offset+40), directory_pdb_B = "%s/pdb/offset_%d_perturbed"%(directory,offset+40), directory_new = directory, parallel = True, n_core = n_core)
    ## Convert PDB to mesh
    convert_pdb_mesh("original", "perturbed", n_sample = n_sample, sm_radius = sm_radius, directory_pdb_A = "%s/pdb/WT_offset_%d"%(directory,offset), directory_pdb_B = "%s/pdb/offset_%d_perturbed"%(directory,offset+40), directory_mesh = "%s/msh_offset_%d/"%(directory,offset), parallel = True, n_core = n_core)

#for n_direction_per_cone in [1,2,4]:
#for cap_radius in [0.30]: #[0.10,0.20,0.40,0.50,0.60,0.70]:
#for n_filtration in [60]: #[20,30,40,50,70,80,90,100,110]:
#for n_cone in [5,10,15,25,30]:
for offset in range(0,50,10):

    print("%d %.2f %d %d %d"%(n_direction_per_cone,cap_radius,n_filtration,n_cone,offset))

    ## Calculate distributed cones of directions for EC calculations
    if not os.path.exists("%s/directions_%d_%d_%.2f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius)):
        directions = generate_equidistributed_cones(n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,hemisphere=False)
        np.savetxt("%s/directions_%d_%d_%.2f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius),directions)
    else:
        directions = np.loadtxt("%s/directions_%d_%d_%.2f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius))
    
    ## EC calculations to convert simplicial meshes to topological summary statistics
    if not os.path.exists("%s/DECT_original_perturbed_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset)):
        X, y, not_vacuum = compute_ec_curve_folder(protA,protB,directions,directory_mesh_A="%s/msh_offset_%d/original_%.1f"%(directory_perturb,offset,sm_radius), directory_mesh_B="%s/msh_offset_%d/perturbed_%.1f"%(directory_perturb,offset,sm_radius),n_sample=n_sample,ec_type=ec_type,n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,n_filtration=n_filtration,sm_radius=sm_radius,parallel=parallel,n_core=n_core) 
        np.savetxt("%s/DECT_original_perturbed_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),X)
        np.savetxt('%s/original_perturbed_label.txt'%directory,y)
        np.savetxt("%s/notvacuum_DECT_original_perturbed_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),not_vacuum)
    else:
        X = np.loadtxt("%s/DECT_original_perturbed_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
        y = np.loadtxt('%s/original_perturbed_label.txt'%directory)
        not_vacuum = np.loadtxt("%s/notvacuum_DECT_original_perturbed_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))

    ## RATE calculation for variable selections from the topological summary statistics
    if not os.path.exists("%s/rates_DECT_original_perturbed_%d_%d_%.2f_%d_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset)):
        kld, rates, delta, eff_samp_size = find_rate_variables_with_other_sampling_methods(X,y,bandwidth=bandwidth,sampling_method=sampling_method,parallel=True,n_core=n_core)
        np.savetxt("%s/rates_DECT_original_perturbed_%d_%d_%.2f_%d_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),rates)
    else:
        rates = np.loadtxt("%s/rates_DECT_original_perturbed_%d_%d_%.2f_%d_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
        
    ## reconstruct the RATE values onto the protein structures for visualization
    ## reconstruct probabilities are stored in "Temperature factor" column in the pdb format
    ## can then be visualized using Chimera or pymol 
    if not os.path.exists("%s/vert_prob_DECT_original_perturbed_%d_%d_%.2f_%d_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset)):
        vert_prob = reconstruct_on_multiple_mesh(protA,protB,directions,rates,not_vacuum,n_sample=n_sample,directory_mesh="%s/msh_offset_%d/original_%.1f/"%(directory_perturb,offset,sm_radius),n_direction_per_cone=n_direction_per_cone,n_filtration=n_filtration,sm_radius=sm_radius)
        np.savetxt("%s/vert_prob_DECT_original_perturbed_%d_%d_%.2f_%d_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),vert_prob)
    else:
        vert_prob = np.loadtxt("%s/vert_prob_DECT_original_perturbed_%d_%d_%.2f_%d_offset_%d.txt"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
    
    if not os.path.exists("%s/vert_prob_DECT_original_perturbed_%d_%d_%.2f_%d_offset_%d.pdb"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset)):
        write_vert_prob_on_pdb(vert_prob,protA=protA,protB=protB,selection=selection,pdb_in_file="WT_R164S_65_213/pdb/WT_offset_%d/%s_frame0.pdb"%(offset,prot),pdb_out_file="%s/vert_prob_DECT_original_perturbed_%d_%d_%.2f_%d_offset_%d.pdb"%(directory,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))

