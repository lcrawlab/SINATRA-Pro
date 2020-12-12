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
ec_type = "DECT"
n_filtration = 25

convert_traj_pdb_aligned(protA, protB, struct_file_A = struct_file_A, traj_file_A = traj_file_A, struct_file_B = struct_file_B, traj_file_B = traj_file_B, nsample = 101, selection = selection)
convert_pdb_mesh(protA, protB, nsample = 101, radius = 4.0)
n_A, n_B, normec_file, notvacuum_file = compute_ec_curve_folder(protA, protB, nsample = nsample, ec_type = "ECT", n_cone = 15, directions_per_cone = 4, cap_radius = 0.1, n_filtration = n_filtration, ball_radius = 1.0, include_faces = True, sm_radius = 4.0)
x = np.loadtxt(normec_file)
N = x.shape[0]
y = np.zeros(N,dtype=int)
y[:n_A].fill(-1)
y[n_A:].fill(1)
kld, rates, delta, eff_samp_size = find_rate_variables_with_other_sampling_methods(x,y,bandwidth=0.01,sampling_method="ESS")
np.savetxt('%s_%s/rates.txt',rates)
print(kld,rates,delta,eff_samp_size)
directions = np.loadtxt('directions_%d_%d_%.1f.txt'%(n_cone,directions_per_cone,cap_radius))

#rates = np.loadtxt(,usecols=1)
not_vacuum = np.loadtxt(not_vacuum)
pip = np.zeros(not_vacuum.size,dtype=float)
j = 0
for i in range(not_vacuum.size):
    if not_vacuum[i]:
        pip[i] = raw[j]
        j += 1

meshProtein = mesh()
meshProtein.read_mesh_file(filename='mesh_noh/%s_4/WT_4_frame0.msh'%protA)
outProb = np.zeros((nframe,meshProtein.vertices.shape[0]),dtype=float)
for i in range(nframe):
    frame = 100*i
    sys.stdout.write('Reconstructing for Frame %d...\r'%frame)
    sys.stdout.flush()
    meshProtein = mesh()
    meshProtein.read_mesh_file(filename='mesh_noh/%s_4/WT_4_frame%d.msh'%(protA,frame))
    outProb[i,:] = reconstruct_by_sorted_threshold_new(meshProtein, directions, pip, n_filtration = nf, n_direction_per_cone = 4, ball_radius = 1.0, by_rank = False)

averageProb = np.average(outProb,axis=0)
np.savetxt('prob_rate_sthrval_vert_%s_%s_%s_%d_15.txt'%(protA,protB,method,nf),averageProb)

prob = np.loadtxt('prob_rate_sthrval_vert_%s_%s_%s_%d_15.txt'%(protA,protB,method,nf))
u = mda.Universe('pdb_noh/%s/WT_frame5000.pdb'%protA)
protein = u.select_atoms('protein and not type H')
u.add_TopologyAttr('tempfactors')
print(np.amin(prob),np.amax(prob))
y = prob*100
print(np.amin(y),np.amax(y))
protein.tempfactors = y
protein.write("%s_%s_%s_%d_15_rate_sthrval.pdb"%(protA,protB,method,nf))
