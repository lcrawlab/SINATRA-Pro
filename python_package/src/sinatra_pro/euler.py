#!/bin/python3

import os, sys
from sinatra_pro.mesh import *
from fast_histogram import histogram1d
import multiprocessing
from joblib import Parallel, delayed

def compute_ec_curve_single(mesh, direction, ball_radius, n_filtration = 25, ec_type = "ECT", include_faces = True):
    """
    Computes the Euler Characteristics (EC) curves in a given direction in discrete filtraion steps for a given mesh 

    `mesh` is the `mesh` class containing vertices, edges, and faces of the mesh.
    
    `direction` is the direction for EC curve to be calculated on.
   
    `ball_radius` is the radius of the bounding ball.
    
    `n_filtration` is the number of sub-level sets for which to compute the EC curve on in a given direction.
    
    `ec_type` is the type of EC transform (ECT), available options: DECT / ECT / SECT. 
    DECT (differential ECT) is the default method used for protein. 
    ECT is the standard ECT and SECT is the smoothe ECT.

    If `included_faces` is set to False, it ignore faces from the EC calculations.
    """
    eulers = np.zeros(n_filtration,dtype=float)
    vertex_function = np.dot(mesh.vertices,direction) 
    radius = np.linspace(-ball_radius,ball_radius,n_filtration)
   
    # filtrating vertices
    V = histogram1d(vertex_function,range=[-ball_radius,ball_radius],bins=(n_filtration-1))
    
    # filtrating edges
    if len(mesh.edges) > 0:
        edge_function = np.amax(vertex_function[mesh.edges],axis=1)
        E = histogram1d(edge_function,range=[-ball_radius,ball_radius],bins=(n_filtration-1))
    else:
        E = 0

    if include_faces and len(mesh.faces) > 0:
        # filtrating faces
        face_function = np.amax(vertex_function[mesh.faces],axis=1)
        F = histogram1d(face_function,range=[-ball_radius,ball_radius],bins=(n_filtration-1))
    else:
        F = 0
    
    eulers[1:] = V - E + F

    if ec_type == "ECT":
        eulers = np.cumsum(eulers)
        return eulers
    elif ec_type == "DECT":
        eulers[1:] = eulers[1:] / (radius[1:]-radius[:-1])
        return eulers
    elif ec_type == "SECT":
        eulers = np.cumsum(eulers)
        eulers -= np.mean(eulers)
        eulers = np.cumsum(eulers)*((radius[-1]-radius[0])/n_filtration)
        return eulers
    else:
        return None

def compute_ec_curve(mesh, directions, n_filtration = 25, ball_radius = 1.0, ec_type = "ECT", first_column_index = False, include_faces = True):
    """Computes the Euler Characteristics (EC) curves in a given direction with single CPU core"""
    eulers = np.zeros((directions.shape[0],n_filtration),dtype=float)
    for i in range(directions.shape[0]):
        eulers[i] = compute_ec_curve_single(mesh,directions[i],ball_radius,n_filtration,ec_type,include_faces)
    radius = np.linspace(-ball_radius,ball_radius,n_filtration)
    return radius, eulers

def compute_ec_curve_parallel(mesh, directions, n_filtration = 25, ball_radius = 1.0, ec_type = "ECT", include_faces = True, n_core = -1):
    """Computes the Euler Characteristics (EC) curves in a given direction with single multiple core"""
    parameter = (ball_radius,n_filtration,ec_type,include_faces)
    if n_core == -1:    
        n_core = multiprocessing.cpu_count()
    processed_list = Parallel(n_jobs=n_core)(delayed(compute_ec_curve_single)(mesh,direction,*parameter) for direction in directions)
    processed_list = np.array(processed_list)
    radius = np.linspace(-ball_radius,ball_radius,n_filtration)
    return radius, processed_list

def compute_ec_curve_folder(protA = "protA", protB = "protB", directions = None, n_sample = 101, ec_type = "ECT", n_filtration = 25, ball_radius = 1.0, include_faces = True, directory_mesh_A = None, directory_mesh_B = None, sm_radius = 4.0, hemisphere=False, parallel = False, n_core = -1, verbose = False):
    """
    Computes the Euler Characteristics (EC) curves for a set of directions for the data set. 
    
    `ec_type` is the type of EC transform (ECT), available options: DECT / ECT / SECT. 
    DECT (differential ECT) is the default method used for protein. 
    ECT is the standard ECT and SECT is the smoothe ECT.

    `directory_mesh_A` and `directory_mesh_B` are the folders that contain the .msh files for meshes for class A and B respectively. 
    If `directory_mesh_A` and `directory_mesh_B` are provided, the function calculates EC curves of all meshes in the two folders.
    
    `directions` is the list of vectors containing all directions for EC curves to be calculated on.
        
    `n_filtration` is the number of sub-level sets for which to compute the EC curve on in a given direction.

    `ball_radius` is the radius for EC calculation, default to be 1.0 for unit sphere, where EC curves are calculated from radius = -1.0 to +1.0.

    If `included_faces` is set to False, it ignore faces from the EC calculations.

    If `parallel` is set to True, the program runs on multiple cores for EC calculations, 
    then `n_core` will be the number of cores used (the program uses all detected cores if `n_core` is not provided`).

    If `verbose` is set to True, the program prints progress in command prompt. 
    """

    if directory_mesh_A == None or directory_mesh_B == None:
        directory = "%s_%s"%(protA,protB)
        directory_mesh_A = "%s_%s/mesh/%s_%.1f"%(protA,protB,protA,sm_radius)
        directory_mesh_B = "%s_%s/mesh/%s_%.1f"%(protA,protB,protB,sm_radius) 
        if not os.path.exists(directory):
            os.mkdir(directory)
        if not os.path.exists(directory_mesh_A):
            os.mkdir(directory_mesh_A)
        if not os.path.exists(directory_mesh_B):
            os.mkdir(directory_mesh_B)

        ecss = []
        for prot, directory_mesh in zip([protA,protB],[directory_mesh_A,directory_mesh_B]):
            ecs = []
            for i_sample in range(n_sample):
                if verbose:
                    sys.stdout.write('Calculating EC for %s Frame %d...\r'%(prot,i_sample))
                    sys.stdout.flush()
                meshProtein = mesh()
                meshProtein.read_mesh_file(filename='%s/%s_frame%d.msh'%(directory_mesh,prot,i_sample))
                if parallel:
                    t, ec = compute_ec_curve_parallel(meshProtein, directions, n_filtration = n_filtration, ball_radius = ball_radius,ec_type = ec_type, include_faces = include_faces, n_core = n_core)
                else:
                    t, ec = compute_ec_curve(meshProtein, directions, n_filtration = n_filtration, ball_radius = ball_radius,ec_type = ec_type, include_faces = include_faces)
                ecs.append(ec.flatten())
            ecs = np.array(ecs)
            ecss.append(ecs)

    else:
        directory = '.'
        ecss = []
        for prot, directory_mesh in zip([protA,protB],[directory_mesh_A,directory_mesh_B]):
            ecs = []
            for filename in os.listdir(directory_mesh):
                if filename.endswith(".msh"):
                    if verbose:
                        sys.stdout.write('Calculating EC for %s %s...\r'%(prot,filename))
                        sys.stdout.flush()
                    meshProtein = mesh()
                    meshProtein.read_mesh_file(filename=directory_mesh + '/' + filename)
                    if parallel:
                        t, ec = compute_ec_curve_parallel(meshProtein, directions, n_filtration = n_filtration, ball_radius = ball_radius,ec_type = ec_type, include_faces = include_faces, n_core = n_core)
                    else:
                        t, ec = compute_ec_curve(meshProtein, directions, n_filtration = n_filtration, ball_radius = ball_radius,ec_type = ec_type, include_faces = include_faces) 
                    ecs.append(ec.flatten())
            ecs = np.array(ecs)
            ecss.append(ecs)
            
    data_A = ecss[0]
    data_B = ecss[1]
    
    vacuum = np.ones(data_A.shape[1],dtype=bool)
    for a in data_A:
        vacuum = np.logical_and(vacuum,a == 0)
    for a in data_B:
        vacuum = np.logical_and(vacuum,a == 0)

    not_vacuum = np.logical_not(vacuum)
    data = np.vstack((data_A,data_B))[:,not_vacuum]

    mean = np.average(data,axis=0)
    std = np.std(data,axis=0)

    data = np.subtract(data,mean)
    data = np.divide(data,std)
    data[np.isnan(data)] = 0
         
    n_A = data_A.shape[0]
    n_B = data_B.shape[0]
    label = np.zeros(n_A+n_B,dtype=int)
    label[:n_A].fill(0)
    label[n_A:].fill(1)

    return data, label, not_vacuum



