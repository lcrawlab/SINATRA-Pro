#!/bin/python3

import os, sys
import numpy as np
import MDAnalysis as mda
from traj_reader import *
from simplices_construction import *
import multiprocessing
from joblib import Parallel, delayed

def calc_radius_write_pdb(directory_original, directory_mesh_A, prot, selection, directory_pdb_B, i_sample):
    msh_file = '%s/%s_frame%d.msh'%(directory_mesh_A,prot,i_sample) 
    original_pdb_file = '%s/%s_frame%d.pdb'%(directory_original,prot,i_sample)
    u = mda.Universe(original_pdb_file).select_atoms('protein')
    original_pos = u.atoms.positions
    r1 = np.amax(np.linalg.norm(u.atoms.positions,axis=1))
    u_new = mda.Universe(original_pdb_file).select_atoms('protein and not %s'%selection)
    u_new.write('%s/%s_frame%d.pdb'%(directory_pdb_B,prot,i_sample))
    r2 = np.amax(np.linalg.norm(u.atoms.positions,axis=1))
    return [r1,r2]

def remove_protein_region(sm_radius = 4.0, selection='resid 163:178', prot = "WT", n_sample = 101, directory_original = 'WT_R164S/pdb/WT/', directory_new = "simulation", directory_pdb_B = None, directory_mesh_A = None, directory_mesh_B = None, parallel = False, n_core = -1):
    
    if parallel:
        if n_core == -1:
            n_core = multiprocessing.cpu_count()   

    if not os.path.exists(directory_new):
        os.mkdir(directory_new)
    directory_mesh = '%s/msh/'%directory_new
    if not os.path.exists(directory_mesh):
        os.mkdir(directory_mesh)
    
    if directory_mesh_A == None:
        directory_mesh_A = '%s/original'%directory_mesh
    if directory_mesh_B == None:
        directory_mesh_B = '%s/perturbed'%directory_mesh
    if not os.path.exists(directory_mesh_A):
        os.mkdir(directory_mesh_A)
    if not os.path.exists(directory_mesh_B):
        os.mkdir(directory_mesh_B)

    if not os.path.exists('%s/pdb'%directory_new):
        os.mkdir('%s/pdb'%directory_new)
    if directory_pdb_B == None:
        directory_pdb_B = '%s/pdb/perturbed/'%directory_new
    if not os.path.exists(directory_pdb_B):
        os.mkdir(directory_pdb_B)

    if parallel:
        processed_list = Parallel(n_jobs=n_core)(delayed(calc_radius_write_pdb)(directory_original, directory_mesh_A, prot, selection, directory_pdb_B, i_sample) for i_sample in range(n_sample))
        r = np.array(processed_list).flatten()
    else:
        r = np.array([])
        for i_sample in range(n_sample):
            output = calc_radius_write_pdb(directory_original, directory_mesh_A, prot, selection, directory_pdb_B, i_sample)
            r = np.append(r,output)
    rmax = np.amax(r)
    
    if parallel:
        Parallel(n_jobs=n_core)(delayed(convert_pdb_mesh_single)(sm_radius=4.0,rmax=rmax,
                                prot=prot,
                                directory_mesh = directory_mesh_A,
                                directory_pdb = directory_original,
                                filename = '%s_frame%d.pdb'%(prot,i_sample),
                                selection='protein') for i_sample in range(n_sample))
        
        Parallel(n_jobs=n_core)(delayed(convert_pdb_mesh_single)(sm_radius=4.0,rmax=rmax,
                                prot=prot,
                                directory_mesh = directory_mesh_B,
                                directory_pdb = directory_pdb_B,
                                filename = '%s_frame%d.pdb'%(prot,i_sample),
                                selection='protein') for i_sample in range(n_sample))
    else:
        
        for i_sample in range(n_sample):
            convert_pdb_mesh_single(sm_radius=4.0,rmax=rmax,
                                prot=prot,
                                directory_mesh = directory_mesh_A,
                                directory_pdb = directory_orignal,
                                filename = '%s_frame%d.pdb'%(prot,i_sample),
                                selection='protein')
            convert_pdb_mesh_single(sm_radius=4.0,rmax=rmax,
                                prot=prot,
                                directory_mesh = directory_mesh_B,
                                directory_pdb = directory_pdb_B,
                                filename = '%s_frame%d.pdb'%(prot,i_sample),
                                selection='protein')

    return

def perturb_protein_region(sm_radius = 4.0, selection='resid 163:178', parallel = False, n_core = -1):
    return

