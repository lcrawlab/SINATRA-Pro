#!/bin/python3

import argparse, os

from sinatra_pro.directions import *
from sinatra_pro.traj_reader import *
from sinatra_pro.euler import *
from sinatra_pro.gp import *
from sinatra_pro.reconstruction import *

##########################################################################

parser = argparse.ArgumentParser(description='SINATRA Pro')

parser.add_argument('-na','--protA', type=str, help='name of protein A for file naming',  default='protA')
parser.add_argument('-nb','--protB', type=str, help='name of protein B for file naming',  default='protB')

parser.add_argument('-sa','--struct_file_A', type=str, help='structure file for protein A (.gro)')
parser.add_argument('-ta','--traj_file_A', type=str, help='trajectory file for protein A (.xtc)')
parser.add_argument('-sb','--struct_file_B', type=str, help='structure file for protein B (.gro)')
parser.add_argument('-tb','--traj_file_B', type=str, help='trajectory file for protein B (.xtc)')

parser.add_argument('-dir','--directory', type=str, help='directory for output files', default='output')

parser.add_argument('-fp' ,'--from_pdb', help='start from sets of PDB structures instead of trajectories', dest='from_pdb', action='store_true')
parser.add_argument('-pa','--pdbpath_A', type=str, help='directory containing PDB structures for protein A')
parser.add_argument('-pb','--pdbpath_B', type=str, help='directory containing PDB structures for protein B')
parser.add_argument('-pr','--pdb_reference', type=str, help='PDB structure for visualization from protein A')

parser.add_argument('-pl','--parallel', dest='parallel', help='use multiple CPU cores for calculations', action='store_true')
parser.add_argument('-nc','--n_core', type=int, help='number of core for parallel computing, default: use all cores', default=-1)

parser.add_argument('-n' ,'--n_sample', type=int, help='number of sample drawn from trajectory, default: 10',default=10)
parser.add_argument('-of','--offset', type=int, help='starting frame for sample drawn from trajectory, default: 0',default=0)
parser.add_argument('-s' ,'--selection', type=str, help='selection for protein, default: all protein', default='protein')
parser.add_argument('-r' ,'--radius', type=float, help='radius for simplicial construction, default: 2.0', default=2.0)
parser.add_argument('-hs','--hemisphere', help='distribute directions over hemisphere instead of whole sphere', dest='hemisphere', action='store_true')

parser.add_argument('-et','--ec_type', type=str, help='type of Euler characteristic measure (DECT/ECT/SECT), default: DECT', default='DECT')
parser.add_argument('-c' ,'--n_cone', type=int, help='number of cone, default: 1', default=1)
parser.add_argument('-d' ,'--n_direction_per_cone', type=int, help='number of direction per cone, default: 1', default=1)
parser.add_argument('-t' ,'--cap_radius', type=float, help='cap radius, default: 0.8', default=0.80)
parser.add_argument('-l' ,'--n_filtration', type=int, help='number of filtration step, default: 20', default=20)

parser.add_argument('-bw','--bandwidth', type=float, help='bandwidth for elliptical slice sampling, default: 0.01',default=0.01)
parser.add_argument('-sm','--sampling_method', type=str, help='sampling method, default: ESS', default='ESS')
parser.add_argument('-nm','--n_mcmc', type=int, help='number of sample from ESS', default=100000)
parser.add_argument('-ll' ,'--logistic_likelihood', help='use logistic likelihood instead of probit likelihood', dest='probit', action='store_false')
parser.add_argument('-lr' ,'--low_rank', help='use low rank matrix approximations to compute the RATE values', dest='low_rank', action='store_true')

parser.add_argument('-v' ,'--verbose', help='verbose', dest='verbose', action='store_true')
parser.add_argument('-no','--name_offset', help='name folder with offset', dest='single', action='store_false')

parser.set_defaults(from_pdb=False,hemisphere=False,probit=True,low_rank=False,parallel=False,verbose=False,single=True)
args = parser.parse_args()

from_pdb = args.from_pdb # if True, start from PDB files

# Name the variants, just for filename purpose
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
hemisphere = args.hemisphere

## Variable selection parameters
bandwidth = args.bandwidth
sampling_method = args.sampling_method
n_mcmc = args.n_mcmc
probit = args.probit
low_rank = args.low_rank
single = args.single
verbose = args.verbose

##########################################################################

## Read trajectory file and output aligned protein structures in pdb format
if not from_pdb:
    convert_traj_pdb_aligned(protA, protB, 
            struct_file_A=struct_file_A, 
            traj_file_A=traj_file_A, 
            struct_file_B=struct_file_B, 
            traj_file_B=traj_file_B, 
            align_frame=offset, 
            n_sample=n_sample, 
            selection=selection, 
            offset=offset, 
            directory=directory,
            single=single, ## single="True" is for single run purpose, "False" for duplicate runs purpose which groups and names file with the frame offset.
            verbose=verbose)

#####################
### IF you already have your own aligned structure, start here
#####################

if not from_pdb:
    directory_pdb_A = "%s/pdb/%s/"%(directory,protA)
    directory_pdb_B = "%s/pdb/%s/"%(directory,protB)
    reference_pdb_file = "%s/%s_frame0.pdb"%(directory_pdb_A,protA) ## which pdb to use for visualization
else:
    # directories for files if start from PDB files
    directory_pdb_A = args.pdbpath_A
    directory_pdb_B = args.pdbpath_B
    reference_pdb_file = args.pdb_reference 

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
        hemisphere=hemisphere)
np.savetxt("%s/directions_%d_%d_%.2f.txt"%(directory,n_cone,n_direction_per_cone,cap_radius),directions)

## EC calculations to convert simplicial meshes to topological summary statistics
X, y, not_vacuum = compute_ec_curve_folder(protA,protB,directions,
        n_sample=n_sample,
        ec_type=ec_type,
        n_filtration=n_filtration,
        sm_radius=sm_radius,
        directory_mesh_A = "%s/msh/%s_%.1f"%(directory,protA,sm_radius), 
        directory_mesh_B = "%s/msh/%s_%.1f"%(directory,protB,sm_radius), 
        parallel=parallel, 
        n_core=n_core, 
        verbose=verbose)

np.savetxt("%s/%s_%s_%s_%.1f_%d_%d_%.2f_%d_norm_all.txt"%(directory,ec_type,protA,protB,sm_radius,n_cone,n_direction_per_cone,cap_radius,n_filtration),X)
np.savetxt("%s/notvacuum_%s_%s_%s_%.1f_%d_%d_%.2f_%d_norm_all.txt"%(directory,ec_type,protA,protB,sm_radius,n_cone,n_direction_per_cone,cap_radius,n_filtration),not_vacuum)
np.savetxt('%s/%s_%s_label_all.txt'%(directory,protA,protB),y)    

## RATE calculation for variable selections from the topological summary statistics
kld, rates, delta, eff_samp_size = calc_rate(X,y,
        bandwidth=bandwidth,
        n_mcmc=n_mcmc,
        low_rank=low_rank,
        parallel=parallel,
        n_core=n_core,
        verbose=verbose)

np.savetxt("%s/rate_%s_%s_%s_%.1f_%d_%d_%.2f_%d.txt"%(directory,ec_type,protA,protB,sm_radius,n_cone,n_direction_per_cone,cap_radius,n_filtration),rates)

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
        directory_mesh="%s/msh/%s_%.1f"%(directory,protA,sm_radius),
        verbose=verbose)

np.savetxt("%s/rate_atom_%s_%s_%s_%.1f_%d_%d_%.2f_%d.txt"%(directory,ec_type,protA,protB,sm_radius,n_cone,n_direction_per_cone,cap_radius,n_filtration),vert_prob)

write_vert_prob_on_pdb(vert_prob,
        protA=protA,
        protB=protB,
        selection=selection, 
        pdb_in_file=reference_pdb_file, 
        pdb_out_file="%s/rate_atom_%s_%s_%s_%.1f_%d_%d_%.2f_%d_all.pdb"%(directory,ec_type,protA,protB,sm_radius,n_cone,n_direction_per_cone,cap_radius,n_filtration))

print("SINATRA Pro calculation completed.")

