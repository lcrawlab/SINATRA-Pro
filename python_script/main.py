#!/bin/python3
<<<<<<< HEAD

=======
>>>>>>> 58eed95a1351010a064e2695526493d115943615
from traj_reader import *
from euler import *
from gp import *
from reconstruction import *
import sys
<<<<<<< HEAD
import argparse

parser = argparse.ArgumentParser(description='SINATRA Pro')

parser.add_argument('-pa','--protA', type=str, help='name of protein A for file naming',  default='protA')
parser.add_argument('-pb','--protB', type=str, help='name of protein B for file naming',  default='protB')

parser.add_argument('-sa','--struct_file_A', type=str, help='structure file for protein A (.gro)')
parser.add_argument('-ta','--traj_file_A', type=str, help='trajectory file for protein A (.xtc)')
parser.add_argument('-sb','--struct_file_B', type=str, help='structure file for protein B (.gro)')
parser.add_argument('-tb','--traj_file_B', type=str, help='trajectory file for protein B (.xtc)')

parser.add_argument('-dir','--directory', type=str, help='directory for output files')

parser.add_argument('-pl','--parallel', dest='parallel', action='store_true')
parser.add_argument('-nc','--n_core', type=int, help='number of core for parallel computing, default: use all cores', default=-1)

parser.add_argument('-n','--n_sample', type=int, help='number of sample drawn from trajectory, default: 10',default=10)
parser.add_argument('-of','--offset', type=int, help='starting frame for sample drawn from trajectory, default: 0',default=0)
parser.add_argument('-s','--selection', type=str, help='selection for protein, default: all protein', default='protein')
parser.add_argument('-r','--radius', type=float, help='radius for simplicial construction, default: 2.0', default=2.0)

parser.add_argument('-et','--ec_type', type=str, help='type of Euler characteristic measure (DECT/ECT/SECT), default: DECT', default='DECT')
parser.add_argument('-c','--n_cone', type=int, help='number of cone, default: 1', default=1)
parser.add_argument('-d','--n_direction_per_cone', type=int, help='number of direction per cone, default: 1', default=1)
parser.add_argument('-t','--cap_radius', type=float, help='cap radius, default: 0.8', default=0.80)
parser.add_argument('-l','--n_filtration', type=int, help='number of filtration step, default: 20', default=20)

parser.add_argument('-bw','--bandwidth', type=float, help='bandwidth for elliptical slice sampling, default: 0.01',default=0.01)
parser.add_argument('-sm','--sampling_method', type=str, help='sampling method, default: ESS', default='ESS')
parser.add_argument('-nm','--n_mcmc', type=int, help='number of sample from ESS', default=100000)

parser.add_argument('-v','--verbose', dest='verbose', action='store_true')

parser.set_defaults(parallel=False,verbose=False)
args = parser.parse_args()

# Protein names
protA = args.protA
protB = args.protB

# Input files
struct_file_A = args.struct_file_A
traj_file_A = args.traj_file_A
struct_file_B = args.struct_file_B
traj_file_B = args.traj_file_B

# Output folder
directory = args.directory

# Parallelization
parallel = args.parallel
n_core = args.n_core

# Sample from trajectory
n_sample = args.n_sample
offset = args.offset
selection = args.selection
sm_radius = args.radius

## EC calculation 
ec_type = args.ec_type
n_cone = args.n_cone
n_direction_per_cone = args.n_direction_per_cone
cap_radius = args.cap_radius
n_filtration = args.n_filtration

## Variable selection parameters
bandwidth = args.bandwidth
sampling_method = args.sampling_method
n_mcmc = args.n_mcmc

verbose = args.verbose

print("%s %s %s %s"%(struct_file_A,traj_file_A,struct_file_B,traj_file_B))
print(directory)
print("%d %d %.1f"%(n_sample,offset,sm_radius))
print("%d %d %.2f %d %d"%(n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))   

## Read trajectory file and output aligned protein structures in pdb format
convert_traj_pdb_aligned(protA, protB, struct_file_A = struct_file_A, traj_file_A = traj_file_A, struct_file_B = struct_file_B, traj_file_B = traj_file_B, align_frame = offset, n_sample = n_sample,selection = selection, offset = offset, directory = directory, verbose = verbose)
## Converted protein structures into simplicial mesehes
convert_pdb_mesh(protA, protB, n_sample = n_sample, sm_radius = sm_radius, directory_pdb_A = "%s/pdb/%s_offset_%d"%(directory,protA,offset), directory_pdb_B = "%s/pdb/%s_offset_%d"%(directory,protB,offset), directory_mesh = "%s/msh_offset_%d/"%(directory,offset), parallel = parallel, n_core = n_core, verbose = verbose)
## Calculate distributed cones of directions for EC calculations
directions = generate_equidistributed_cones(n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,hemisphere=False)
np.savetxt("%s/directions_%d_%d_%.2f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius),directions)

## EC calculations to convert simplicial meshes to topological summary statistics
X, y, not_vacuum = compute_ec_curve_folder(protA,protB,directions,n_sample=n_sample,ec_type=ec_type,n_cone=n_cone,n_direction_per_cone=n_direction_per_cone,cap_radius=cap_radius,n_filtration=n_filtration,sm_radius=sm_radius,directory_mesh_A = "%s/msh_offset_%d/%s_%.1f"%(directory,offset,protA,sm_radius), directory_mesh_B = "%s/msh_offset_%d/%s_%.1f"%(directory,offset,protB,sm_radius), parallel=parallel,n_core=n_core,verbose=verbose)
np.savetxt("%s/%s_%s_%s_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),X)
np.savetxt("%s/notvacuum_%s_%s_%s_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),not_vacuum)
np.savetxt('%s/%s_%s_label.txt'%(directory,protA,protB),y)    
#X = np.loadtxt("%s/%s_%s_%s_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
#not_vacuum = np.loadtxt("%s/notvacuum_%s_%s_%s_%d_%d_%.2f_%d_norm_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
#y = np.loadtxt('%s/%s_%s_label.txt'%(directory,protA,protB))

## RATE calculation for variable selections from the topological summary statistics
kld, rates, delta, eff_samp_size = find_rate_variables_with_other_sampling_methods(X,y,bandwidth=bandwidth,sampling_method=sampling_method,n_mcmc=n_mcmc,parallel=parallel,n_core=n_core,verbose=verbose)
np.savetxt("%s/rates_%s_%s_%s_%d_%d_%.2f_%d_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),rates)
#rates = np.loadtxt("%s/rates_%s_%s_%s_%d_%d_%.2f_%d_offset_%d.txt"%(directory,ec_type,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))
=======

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
>>>>>>> 58eed95a1351010a064e2695526493d115943615

## reconstruct the RATE values onto the protein structures for visualization
## reconstruct probabilities are stored in "Temperature factor" column in the pdb format
## can then be visualized using Chimera or Pymol  
<<<<<<< HEAD
vert_prob = reconstruct_on_multiple_mesh(protA,protB,directions,rates,not_vacuum,n_sample=n_sample,n_direction_per_cone=n_direction_per_cone,n_filtration=n_filtration,sm_radius=sm_radius,directory_mesh="%s/msh_offset_%d/%s_%.1f"%(directory,offset,protA,sm_radius),verbose=verbose)
np.savetxt("%s/vert_prob_DECT_%s_%s_%d_%d_%.2f_%d_offset_%d.txt"%(directory,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset),vert_prob)
write_vert_prob_on_pdb(vert_prob,protA=protA,protB=protB,selection=selection,pdb_in_file = "%s/pdb/%s_offset_%d/%s_frame0.pdb"%(directory,protA,offset,protA),pdb_out_file="%s/vert_prob_DECT_%s_%s_%d_%d_%.2f_%d_offset_%d.pdb"%(directory,protA,protB,n_cone,n_direction_per_cone,cap_radius,n_filtration,offset))


=======
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
>>>>>>> 58eed95a1351010a064e2695526493d115943615
