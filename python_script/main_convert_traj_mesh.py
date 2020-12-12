#!/bin/python3

from traj_reader import *

protA = "WT"
protB = "R164S"
nsample = 101

struct_file_A = '%s/md_0_1.gro'%protA
traj_file_A = '%s/md_0_1_noPBC.xtc'%protA

struct_file_B = '%s/md_0_1.gro'%protB
traj_file_B = '%s/md_0_1_noPBC.xtc'%protB

selection = 'resid 65:213'

convert_traj_pdb_aligned(protA, protB, struct_file_A = struct_file_A, traj_file_A = traj_file_A, struct_file_B = struct_file_B, traj_file_B = traj_file_B, nsample = 101, selection = selection):
convert_pdb_mesh(protA, protB, nsample = 101, radius = 4.0)

