#!/bin/python3

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import numpy as np, os, sys
from simplices_construction import *

####
# Convert MD simulation trajectory to aligned protein structures in PDB format
# 
# protA = name of protein A
# protB = name of protein B
# struct_file = protein structure file e.g. .pdb, .gro, .tpr, .pqr
# traj_file = MD trajectory file e.g. .xtc, .dcd, .trj
# n_sample = # of sample structure drawn from trajecotry at even time interval
# directory = directory to store all the files
####

def convert_traj_pdb_aligned(protA, protB, struct_file_A, traj_file_A, struct_file_B, traj_file_B, align_frame = 0, n_sample = 100, selection = None, directory = None, offset = 0):
    if directory == None:
        directory = "%s_%s"%(protA,protB)

    if not os.path.exists(directory):
        os.mkdir(directory)
    
    if not os.path.exists("%s/pdb"%directory):
        os.mkdir("%s/pdb"%directory)

    refu = mda.Universe(struct_file_A,traj_file_A)
    refu.trajectory[align_frame]
    if selection == None:
        refuCA = refu.select_atoms('name CA')
    else:
        refuCA = refu.select_atoms('name CA and %s'%selection)
    com_refuCA = refuCA.center_of_mass()
    refu0 = refuCA.positions - com_refuCA

    for prot, struct_file, traj_file in zip([protA,protB],[struct_file_A,struct_file_B],[traj_file_A,traj_file_B]):
        directory_pdb = "%s/pdb/%s_offset_%d"%(directory,prot,offset)
        if not os.path.exists(directory_pdb):
            os.mkdir(directory_pdb)

        u = mda.Universe(struct_file,traj_file)
        u.trajectory[align_frame]
        if selection == None:
            refCA = refu.select_atoms('name CA')
        else:
            refCA = u.select_atoms('name CA and %s'%selection)

        com_refCA = refCA.center_of_mass()
        ref0 = refCA.positions - com_refCA

        R, rmsdval = align.rotation_matrix(ref0, refu0)
        refCA.translate(-refCA.center_of_mass())
        refCA.rotate(R)
        ref0 = refCA.positions

        if selection == None:
            CA = u.select_atoms('name CA')
            atoms = u.select_atoms('protein')
        else:
            CA = u.select_atoms('name CA and %s'%selection)
            atoms = u.select_atoms('protein and %s'%selection)

        rmsds = []
        t = []
        nframe = len(u.trajectory)
        nskip =  int(nframe/n_sample)
        
        frame = 0
        i_sample = 0
        for ts in u.trajectory:
            if (frame-offset) % nskip == 0:
                sys.stdout.write("Writing pdb files for %s, t = %.1f\r"%(prot,ts.time))
                sys.stdout.flush()
                traj0 = CA.positions - CA.center_of_mass()
                R, rmsdval = align.rotation_matrix(traj0, ref0)
                atoms.translate(-CA.center_of_mass())
                atoms.rotate(R)
                atoms.write('%s/%s_frame%d.pdb'%(directory_pdb,prot,i_sample))
                i_sample += 1
                if i_sample == n_sample:
                    break
            frame += 1
        sys.stdout.write("\n")
    return

####
# Convert the aligned protein structures in PDB format (e.g. from "convert_traj_pdb_aligned") to simplicial meshes
# 
# protA = name of protein A
# protB = name of protein B 
# n_sample = # of sample structure drawn from trajecotry at even time interval
# radius = cutoff radius for constructing simplicial meshes
# directory_pdb_A = directory for the input pdb files, default = protA_protB/pdb/protA if not specified
# directory_pdb_B = directory for the input pdb files, default = protA_protB/pdb/protB if not specified
# directory_mesh = directory for the input pdb files, default = protA_protB/mesh/ if not specified
####

def calc_radius_pdb(selection='protein',directory=None,prot=None,i_sample=None,directory_pdb=None,filename=None):
    if directory != None and prot != None and i_sample != None:
        pdb_file = '%s/pdb/%s/%s_frame%d.pdb'%(directory,prot,prot,i_sample)
    if directory_pdb != None:
        if not filename.endswith(".pdb"):
            return None
        else:
            pdb_file = directory_pdb + '/' + filename
    trajwhole = mda.Universe(pdb_file).select_atoms(selection)
    comp = ComplexFiltration()
    comp.vertices = trajwhole.positions
    return comp.calc_radius()

def convert_pdb_mesh_single(sm_radius, rmax, directory = None, prot = None , i_sample = None, directory_mesh = None, directory_pdb = None, filename = None, selection='protein'):
    if directory != None and prot != None and i_sample != None and directory_mesh != None:
        pdb_file = '%s/pdb/%s/%s_frame%d.pdb'%(directory,prot,prot,i_sample)
        msh_file = '%s/%s_frame%d.msh'%(directory_mesh,prot,i_sample)
        sys.stdout.write('Constructing topology for %s for Frame %d...\r'%(prot,i_sample))
        sys.stdout.flush()
    if directory_mesh != None and directory_pdb != None and filename != None and prot != None:
        pdb_file = directory_pdb + '/' + filename
        msh_file = '%s/%s.msh'%(directory_mesh,filename[:-4])
        sys.stdout.write('Constructing topology for %s for %s...\r'%(prot,filename))
        sys.stdout.flush()
    trajwhole = mda.Universe(pdb_file).select_atoms(selection)
    comp = ComplexFiltration()
    comp.vertices = trajwhole.positions
    comp.calc_distance_matrix()
    edges, distances = comp.get_edge_list(radius=sm_radius)
    faces = comp.edge_to_face_list(edges=edges)
    comp.vertices /= rmax
    comp.write_mesh_file(edges=edges,faces=faces,filename=msh_file)
    return


def convert_pdb_mesh(protA, protB, n_sample = 101, sm_radius = 4.0, directory_pdb_A = None, directory_pdb_B = None, directory_mesh = None, parallel = False, n_core = -1):
    
    if parallel:
        import multiprocessing
        from joblib import Parallel, delayed
        if n_core == -1:    
            n_core = multiprocessing.cpu_count()

    if directory_pdb_A == None or directory_pdb_B == None:
        directory = "%s_%s"%(protA,protB)
        if not os.path.exists(directory):
            os.mkdir(directory)
        directory_mesh = "%s/mesh"%directory
        if not os.path.exists(directory_mesh):
            os.mkdir(directory_mesh)
    else:
        if directory_mesh == None:
            directory = '.' 
            directory_mesh = "mesh"
            if not os.path.exists("mesh"):
                os.mkdir("mesh")
        else:
            if not os.path.exists(directory_mesh):
                os.mkdir(directory_mesh)


    r = np.array([])
    sys.stdout.write('Calculating sphere radius...\n')
    if directory_pdb_A == None or directory_pdb_B == None:
        for prot in [protA, protB]:
            if parallel:
                r_prot = Parallel(n_jobs=n_core)(delayed(calc_radius_pdb)(selection='protein',directory=directory,prot=prot,i_sample=i_sample) for i_sample in range(n_sample))
                r = np.append(r,r_prot)
            else:
                for i_sample in range(n_sample):
                    r_prot = calc_radius_pdb(selection='protein',directory=directory,prot=prot,i_sample=i_sample)
                    r = np.append(r,r_prot)
    else:
        for directory_pdb in [directory_pdb_A,directory_pdb_B]:
            if parallel:
                r_prot = Parallel(n_jobs=n_core)(delayed(calc_radius_pdb)(selection='protein',directory_pdb=directory_pdb,filename=filename) for filename in os.listdir(directory_pdb))
                r = np.append(r,r_prot)
            else:
                for filename in os.listdir(directory_pdb):
                    r_prot = calc_radius_pdb(selection='protein',directory_pdb=directory_pdb,filename=filename)
                    if r_prot != None:
                        r = np.append(r,r_prot)

    rmax = np.amax(r)
    sys.stdout.write('Rmax = %.3f\n'%rmax)
    
    if directory_pdb_A == None or directory_pdb_B == None:
        for prot in [protA, protB]:
            directory_mesh = "%s/mesh/%s_%.1f/"%(directory,prot,sm_radius)
            if not os.path.exists(directory_mesh):
                os.mkdir(directory_mesh)
            if parallel:
                tmp = Parallel(n_jobs=n_core)(delayed(convert_pdb_mesh_single)(sm_radius=sm_radius,rmax=rmax,directory=directory,prot=prot,i_sample=i_sample,directory_mesh=directory_mesh,selection='protein') for i_sample in range(n_sample)) 
            else:
                for i_sample in range(n_sample): 
                    convert_pdb_mesh_single(sm_radius=sm_radius,rmax=rmax,directory=directory,prot=prot,i_sample=i_sample,directory_mesh=directory_mesh,selection='protein')
            sys.stdout.write('\n')
    else:
        for prot, directory_pdb in zip([protA,protB],[directory_pdb_A,directory_pdb_B]):
            directory_mesh_prot = '%s/%s_%.1f'%(directory_mesh,prot,sm_radius)
            if not os.path.exists(directory_mesh_prot):
                os.mkdir(directory_mesh_prot)
            if parallel:
                tmp = Parallel(n_jobs=n_core)(delayed(convert_pdb_mesh_single)(sm_radius=sm_radius,rmax=rmax,directory_pdb=directory_pdb,filename=filename,directory_mesh=directory_mesh_prot,prot=prot,selection='protein') for filename in os.listdir(directory_pdb))
            else:
                for filename in os.listdir(directory_pdb):
                     convert_pdb_mesh_single(sm_radius=sm_radius,rmax=rmax,directory_pdb=directory_pdb,filename=filename,directory_mesh=directory_mesh,prot=prot,selection='protein')
            sys.stdout.write('\n') 
    return

def convert_traj_pdb_single_offset(prot, struct_file, traj_file, align_frame = 0, n_sample = 100, offset = 0, selection = None, directory = None):
    if directory == None:
        directory = "%s"%prot

    if not os.path.exists(directory):
        os.mkdir(directory)
    
    if not os.path.exists("%s/pdb"%directory):
        os.mkdir("%s/pdb"%directory)

    refu = mda.Universe(struct_file,traj_file)
    refu.trajectory[align_frame]
    if selection == None:
        refuCA = refu.select_atoms('name CA')
    else:
        refuCA = refu.select_atoms('name CA and %s'%selection)
    com_refuCA = refuCA.center_of_mass()
    refu0 = refuCA.positions - com_refuCA

    directory_pdb = "%s/pdb/offset_%d"%(directory,offset)
    if not os.path.exists(directory_pdb):
        os.mkdir(directory_pdb)

    u = mda.Universe(struct_file,traj_file)
    u.trajectory[align_frame]
    if selection == None:
        refuCA = refu.select_atoms('name CA')
    else:
        refCA = u.select_atoms('name CA and %s'%selection)

    ref0 = refCA.positions

    CA = u.select_atoms('name CA and %s'%selection)
    if selection == None:
        atoms = u.select_atoms('protein')
    else:
        atoms = u.select_atoms('protein and %s'%selection)

    rmsds = []
    t = []
    nframe = len(u.trajectory)
    nskip =  int(nframe/n_sample)
    
    frame = 0
    i_sample = 0
    for ts in u.trajectory:
        if (frame - offset) % nskip == 0:
            sys.stdout.write("Writing pdb files for %s, offset = %d, t = %.1f\r"%(prot,offset,ts.time))
            sys.stdout.flush()
            traj0 = CA.positions - CA.center_of_mass()
            R, rmsdval = align.rotation_matrix(traj0, ref0)
            atoms.translate(-CA.center_of_mass())
            atoms.rotate(R)
            atoms.write('%s/%s_frame%d.pdb'%(directory_pdb,prot,i_sample))
            i_sample += 1
            if i_sample == n_sample:
                break
        frame += 1
    sys.stdout.write("\n")
    return

