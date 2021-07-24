#!/bin/python3

import os, sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
from MDAnalysis.core.groups import AtomGroup
from sinatra_pro.mesh import *

def convert_traj_pdb_aligned(protA, protB, struct_file_A, traj_file_A, struct_file_B, traj_file_B, align_frame = 0, n_sample = 100, selection = None, directory = None, offset = 0, align_sequence = False, single = False, verbose = False):
    """
    Convert MD simulation trajectory to aligned protein structures in PDB format
    
    Protein structure alignment is done by finding the rotational matrix that minimize the root mean square distance (RMSD) between two comparing coordinate sets.

    `protA` and `protB` are the name of protein A and B respectively for file naming.

    `struct_file_A` and `struct_file_B` are the the protein structure file (e.g. .pdb, .gro, .tpr, .pqr.) for protein A and B respectively.

    `traj_file_A` and `traj_file_B` are the Molecular Dynamics (MD) simulation trajectory file (e.g. .xtc, .dcd, .trj) for protein A and B respectively.

    `n_sample` is the number of sample structure drawn from trajecotry at even time interval.
    
    `directory` is the directory to store the output files.
    
    `offset` is the frame number to start drawing samples from the trajectory.

    If `align_sequence` is set to True, the program aligns the sequence using the Needleman-Wunsch algorithm implemented by MDAnalysis.

    If `single` is set to True, the program just names the folders containing the output files without the `offset` in the folder name.

    If `verbose` is set to True, the program prints progress in command prompt.
    """
   
    if directory == None:
        directory = "%s_%s"%(protA,protB)

    if not os.path.exists(directory):
        os.mkdir(directory)
    
    if not os.path.exists("%s/pdb"%directory):
        os.mkdir("%s/pdb"%directory)
    
    seqselA = []
    seqselB = []

    if selection == None:
        selection = 'protein'

    if align_sequence:
        u_A = mda.Universe(struct_file_A).select_atoms(selection)
        u_B = mda.Universe(struct_file_B).select_atoms(selection)
        seq_align = align.sequence_alignment(u_A,u_B)
        seqA = seq_align[1]
        seqB = seq_align[0]
        nres = seq_align[4]
        for i in range(nres):
            if seqA[i] != '-' and seqB[i] != '-':
                seqselA.append(True)
                seqselB.append(True)
            elif seqA[i] == '-' and seqB[i] != '-':
                seqselB.append(False)
            elif seqB[i] == '-' and seqA[i] != '-':
                seqselA.append(False)
    refu = mda.Universe(struct_file_A,traj_file_A)
    refu.trajectory[align_frame]

    if align_sequence:
        protein = refu.select_atoms(selection)
        atoms_A = AtomGroup([],refu)
        for i, a in enumerate(protein.residues):
            if seqselA[i]:
                atoms_A = atoms_A + a.atoms
        groundref = atoms_A
        groundrefCA = atoms_A.select_atoms('name CA')
    else:
        groundref = refu.select_atoms(selection)
        groundrefCA = refu.select_atoms('name CA and %s'%selection)

    groundref.translate(-groundrefCA.center_of_mass())
    
    for prot, seqsel, struct_file, traj_file in zip([protA,protB],[seqselA,seqselB],[struct_file_A,struct_file_B],[traj_file_A,traj_file_B]):
            
        if not single:    
            directory_pdb = "%s/pdb/%s_offset_%d"%(directory,prot,offset)
        else:
            directory_pdb = "%s/pdb/%s"%(directory,prot)

        if not os.path.exists(directory_pdb):
            os.mkdir(directory_pdb)

        u = mda.Universe(struct_file,traj_file)
        u.trajectory[align_frame]
        
        rmsds = []
        t = []
        n_frame = len(u.trajectory)
        
        if n_sample > n_frame:
            print("n_sample > number of frames in trajectory files!")
            exit()
        
        nskip =  int(n_frame/n_sample)
         
        frame = 0
        i_sample = 0
        for ts in u.trajectory:
            if (frame-offset) % nskip == 0:
                if verbose:
                    sys.stdout.write("Writing pdb files for %s, t = %.1f\r"%(prot,ts.time))
                    sys.stdout.flush()
                protein = u.select_atoms(selection)
                if align_sequence:
                    mobile = AtomGroup([],u)
                    for i, a in enumerate(protein.residues):
                        if seqsel[i]:
                            mobile = mobile + a.atoms
                else:
                    mobile = protein.atoms
                align.alignto(mobile,groundref,select="name CA",weights="mass")
                mobile.atoms.write('%s/%s_frame%d.pdb'%(directory_pdb,prot,i_sample))
                i_sample += 1
                if i_sample == n_sample:
                    break
            frame += 1
        if verbose:
            sys.stdout.write("\n") 
    return

def calc_radius_pdb(selection='protein',directory=None,prot=None,i_sample=None,directory_pdb=None,filename=None):
    """Calculate radius of a PDB structure i.e. distance of the atom furthest away from origin."""
    if directory != None and prot != None and i_sample != None:
        pdb_file = '%s/pdb/%s/%s_frame%d.pdb'%(directory,prot,prot,i_sample)
    if directory_pdb != None:
        if not filename.endswith(".pdb"):
            return None
        else:
            pdb_file = directory_pdb + '/' + filename
    protein = mda.Universe(pdb_file).select_atoms(selection)
    if len(protein) == 0:
        print("PDB file %s selection %s is empty!"%(pdb_file,selection))
        exit()
    meshA = mesh()
    meshA.vertices = protein.positions
    return meshA.calc_radius()

def convert_pdb_mesh_single(sm_radius, rmax, directory = None, prot = None , i_sample = None, directory_mesh = None, directory_pdb = None, filename = None, selection='protein', verbose = False):
    """Convert PDB structure to mesh by simplical construction."""
    if directory != None and prot != None and i_sample != None and directory_mesh != None:
        pdb_file = '%s/pdb/%s/%s_frame%d.pdb'%(directory,prot,prot,i_sample)
        msh_file = '%s/%s_frame%d.msh'%(directory_mesh,prot,i_sample)
        if verbose:
            sys.stdout.write('Constructing topology for %s for Frame %d...\r'%(prot,i_sample))
            sys.stdout.flush()
    if directory_mesh != None and directory_pdb != None and filename != None and prot != None:
        pdb_file = directory_pdb + '/' + filename
        msh_file = '%s/%s.msh'%(directory_mesh,filename[:-4])
        if verbose:
            sys.stdout.write('Constructing topology for %s for %s...\r'%(prot,filename))
            sys.stdout.flush()
    u = mda.Universe(pdb_file)
    protein = u.select_atoms(selection)
    meshA = mesh()
    meshA.vertices = protein.positions
    meshA.convert_vertices_to_mesh(sm_radius=sm_radius,msh_file=msh_file,rmax=rmax)
    return

def convert_pdb_mesh(protA = "protA", protB = "protB", n_sample = 101, sm_radius = 4.0, directory_pdb_A = None, directory_pdb_B = None, directory_mesh = None, parallel = False, n_core = -1, verbose = False):
    """
    Convert the aligned protein structures in PDB format (e.g. from "convert_traj_pdb_aligned") to simplicial meshes
    
    `radius` is the cutoff radius for constructing simplicial meshes.
    
    `directory_pdb_A` is the directory for the input pdb files, default = protA_protB/pdb/protA if not specified.
    
    `directory_pdb_B` is the directory for the input pdb files, default = protA_protB/pdb/protB if not specified.
    
    `directory_mesh` is the directory for the input pdb files, default = protA_protB/mesh/ if not specified.
    """

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
    if verbose:
        sys.stdout.write('Calculating sphere radius...\n')
    if directory_pdb_A == None or directory_pdb_B == None:
        for prot in [protA, protB]:
            if parallel:
                r_prots = Parallel(n_jobs=n_core)(delayed(calc_radius_pdb)(selection='protein',directory=directory,prot=prot,i_sample=i_sample) for i_sample in range(n_sample))
                r = np.append(r,r_prots)
            else:
                for i_sample in range(n_sample):
                    r_prots = calc_radius_pdb(selection='protein',directory=directory,prot=prot,i_sample=i_sample)
                    r = np.append(r,r_prots)
    else:
        for directory_pdb in [directory_pdb_A,directory_pdb_B]:
            if parallel:
                r_prots = Parallel(n_jobs=n_core)(delayed(calc_radius_pdb)(selection='protein',directory_pdb=directory_pdb,filename=filename) for filename in os.listdir(directory_pdb))
                r = np.append(r,r_prots)
            else:
                r_prots = []
                for filename in os.listdir(directory_pdb):
                    r_prot = calc_radius_pdb(selection='protein',directory_pdb=directory_pdb,filename=filename)
                    r_prots.append(r_prot)
                r = np.append(r,r_prots)
            if len(r_prots) == 0:
                print("Folder %s is emply or contain no PDB files!"%directory_pdb)
                exit()

    rmax = np.amax(r)
    if verbose:
        sys.stdout.write('Rmax = %.3f\n'%rmax)
    
    if directory_pdb_A == None or directory_pdb_B == None:
        for prot in [protA, protB]:
            directory_mesh = "%s/mesh/%s_%.1f/"%(directory,prot,sm_radius)
            if not os.path.exists(directory_mesh):
                os.mkdir(directory_mesh)
            if parallel:
                tmp = Parallel(n_jobs=n_core)(delayed(convert_pdb_mesh_single)(sm_radius=sm_radius,rmax=rmax,directory=directory,prot=prot,i_sample=i_sample,directory_mesh=directory_mesh,selection='protein',verbose=verbose) for i_sample in range(n_sample)) 
            else:
                for i_sample in range(n_sample): 
                    convert_pdb_mesh_single(sm_radius=sm_radius,rmax=rmax,directory=directory,prot=prot,i_sample=i_sample,directory_mesh=directory_mesh,selection='protein',verbose=verbose)
    else:
        for prot, directory_pdb in zip([protA,protB],[directory_pdb_A,directory_pdb_B]):
            directory_mesh_prot = '%s/%s_%.1f/'%(directory_mesh,prot,sm_radius)
            if not os.path.exists(directory_mesh_prot):
                os.mkdir(directory_mesh_prot)
            if parallel:
                tmp = Parallel(n_jobs=n_core)(delayed(convert_pdb_mesh_single)(sm_radius=sm_radius,rmax=rmax,directory_pdb=directory_pdb,filename=filename,directory_mesh=directory_mesh_prot,prot=prot,selection='protein',verbose=verbose) for filename in os.listdir(directory_pdb))
            else:
                for filename in os.listdir(directory_pdb):
                    convert_pdb_mesh_single(sm_radius=sm_radius,rmax=rmax,directory_pdb=directory_pdb,filename=filename,directory_mesh=directory_mesh_prot,prot=prot,selection='protein',verbose=verbose)
    if verbose:
        sys.stdout.write('\n')

    return

