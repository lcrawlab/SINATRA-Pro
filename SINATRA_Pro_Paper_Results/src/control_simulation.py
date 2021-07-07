#!/bin/python3

import os, sys
import numpy as np
import MDAnalysis as mda
from traj_reader import *
from simplices_construction import *
import multiprocessing
from joblib import Parallel, delayed


def calc_radius_write_pdb_perturb_both(directory_original, prot, selection, directory_pdb_B, i_sample, atom_noise = 0.1, dr = 1.0, perturb_noise = 0.5, dr_vector = [0.5,0.5,0.5]):
    original_pdb_file = '%s/%s_frame%d.pdb'%(directory_original,prot,i_sample)
    u = mda.Universe(original_pdb_file).select_atoms('protein')
    original_pos = u.atoms.positions
    r1 = np.amax(np.linalg.norm(u.atoms.positions,axis=1))
    u_new = mda.Universe(original_pdb_file)
    region = u_new.select_atoms(selection).ix
    if atom_noise > 0.0:
        u_new.atoms.positions[:,0] += np.random.normal(loc=0.0,scale=atom_noise,size=len(u.atoms))
        u_new.atoms.positions[:,1] += np.random.normal(loc=0.0,scale=atom_noise,size=len(u.atoms))
        u_new.atoms.positions[:,2] += np.random.normal(loc=0.0,scale=atom_noise,size=len(u.atoms))
    if perturb_noise > 0.0:
        perturb_noises = np.random.normal(loc=0.0,scale=perturb_noise,size=len(region))
        j = 0
    for i in range(len(u.atoms)):
        if u_new.atoms[i].ix in region:
            r = np.linalg.norm(u_new.atoms[i].position)
            if dr > 0.0:
                if perturb_noise > 0.0:
                    u_new.atoms[i].position *= (r+dr+perturb_noises[j])/r
                    j += 1
                else:
                    u_new.atoms[i].position *= (r+dr)/r
            u_new.atoms[i].position += dr_vector
    u_new.atoms.write('%s/%s_frame%d.pdb'%(directory_pdb_B,prot,i_sample))
    r2 = np.amax(np.linalg.norm(u.atoms.positions,axis=1))
    return [r1,r2]


def perturb_protein_region_both(sm_radius = 4.0, selection='resid 163:178', atom_noise = 0.5, dr = 0.5, perturb_noise = 0.5, dr_vector = [0.5,0.5,0.5], prot = "WT", n_sample = 101, directory_original = 'WT_R164S/pdb/WT/', directory_new = "simulation", directory_pdb_B = None, parallel = False, n_core = -1):   
 
    if parallel:
        if n_core == -1:
            n_core = multiprocessing.cpu_count()   

    if not os.path.exists(directory_new):
        os.mkdir(directory_new)    
    if not os.path.exists('%s/pdb'%directory_new):
        os.mkdir('%s/pdb'%directory_new)
    if directory_pdb_B == None:
        directory_pdb_B = '%s/pdb/perturbed/'%directory_new
    if not os.path.exists(directory_pdb_B):
        os.mkdir(directory_pdb_B)

    if parallel:
        processed_list = Parallel(n_jobs=n_core)(delayed(calc_radius_write_pdb_perturb_both)(directory_original, prot, selection, directory_pdb_B, i_sample, atom_noise, dr, perturb_noise, dr_vector) for i_sample in range(n_sample))
        r = np.array(processed_list).flatten()
    else:
        r = np.array([])
        for i_sample in range(n_sample):
            output = calc_radius_write_pdb_perturb_new(directory_original, prot, selection, directory_pdb_B, i_sample, atom_noise, dr, perturb_noise, dr_vector)
            r = np.append(r,output)
    rmax = np.amax(r)
    
    return rmax

def calc_radius_write_pdb_perturb_sphere(directory_original, prot, selection, directory_pdb_B, i_sample, r = 1.0, noise = 0.0):
    original_pdb_file = '%s/%s_frame%d.pdb'%(directory_original,prot,i_sample)
    u = mda.Universe(original_pdb_file).select_atoms('protein')
    original_pos = u.atoms.positions
    r1 = np.amax(np.linalg.norm(u.atoms.positions,axis=1))
    u_new = mda.Universe(original_pdb_file)
    region = u_new.select_atoms(selection).ix
    
    x = np.random.normal(loc=0.0,scale=1.0,size=len(region))
    y = np.random.normal(loc=0.0,scale=1.0,size=len(region))
    z = np.random.normal(loc=0.0,scale=1.0,size=len(region))
    v = np.vstack((x,y,z))
    if noise > 0.0:
         r = np.random.normal(loc=r,scale=noise,size=len(region))
    v = v/np.linalg.norm(v,axis=0)*r
    v = v.T
    k = 0
    for i in range(len(u.atoms)):
        if u_new.atoms[i].ix in region:
            u_new.atoms[i].position += v[k]
            k += 1
    u_new.atoms.write('%s/%s_frame%d.pdb'%(directory_pdb_B,prot,i_sample))
    r2 = np.amax(np.linalg.norm(u.atoms.positions,axis=1))
    return [r1,r2]


def perturb_protein_region_sphere(sm_radius = 4.0, selection='resid 163:178', r = 1.0, noise = 0.0, prot = "WT", n_sample = 101, directory_original = 'WT_R164S/pdb/WT/', directory_new = "simulation", directory_pdb_B = None, parallel = False, n_core = -1):   
 
    if parallel:
        if n_core == -1:
            n_core = multiprocessing.cpu_count()
    if not os.path.exists(directory_new):
        os.mkdir(directory_new)
    if not os.path.exists('%s/pdb'%directory_new):
        os.mkdir('%s/pdb'%directory_new)
    if directory_pdb_B == None:
        directory_pdb_B = '%s/pdb/perturbed/'%directory_new
    if not os.path.exists(directory_pdb_B):
        os.mkdir(directory_pdb_B)

    if parallel:
        processed_list = Parallel(n_jobs=n_core)(delayed(calc_radius_write_pdb_perturb_sphere)(directory_original, prot, selection, directory_pdb_B, i_sample, r, noise) for i_sample in range(n_sample))
        r = np.array(processed_list).flatten()
    else:
        r = np.array([])
        for i_sample in range(n_sample):
            output = calc_radius_write_pdb_perturb_new(directory_original, prot, selection, directory_pdb_B, i_sample, r, noise)
            r = np.append(r,output)
    rmax = np.amax(r)
    return rmax

