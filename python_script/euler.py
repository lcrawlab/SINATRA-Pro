#!/bin/python3

from mesh import *
from directions import *
import numpy as np
import multiprocessing
from joblib import Parallel, delayed

#' Computes the EC curve on a ball
#'
#' @export
#' @description Computes the EC curve of a mesh on a bounding ball of specified radius $radius$. The filtration steps for which the
#' (S/D) EC curve are computed is done relative to the bounding ball. For filtrations steps that are computed relative to the shape, use
#' \code{compute_discrete_ec_curve}. For comparisons with multiple shapes, this is the recommended function.
#'
#' @param complex (list) : The list containing metadata about the mesh
#' @param directions (array)
#' @param n_filtration (int): The number of sub-level sets for which to compute the EC curve on in a given direction.
#' @param ball_radius (float): The radius of the bounding ball.
#'
#' @return radius
#' @return ec_curve (nxp matrix): A matrix containing the EC curve of the given simplicial complex with the first index as the projections for which the EC was computed on.

def compute_ec_curve_single(mesh, direction, radius, n_filtration = 25, ec_type = "ECT", include_faces = True):
    eulers = np.zeros(n_filtration,dtype=float)
    vertex_function = np.dot(mesh.vertices,direction) 
    # filtrating vertices
    vertex_hist, bin_edges = np.histogram(vertex_function,bins=radius)
    V = np.cumsum(vertex_hist)
    # filtrating edges
    edge_function = np.amax(vertex_function[mesh.edges],axis=1)
    edge_hist, bin_edges = np.histogram(edge_function,bins=radius)
    E = np.cumsum(edge_hist)
    if include_faces:
        # filtrating faces
        face_function = [np.amax(vertex_function[face]) for face in mesh.faces]
        face_hist, bin_edges = np.histogram(face_function,bins=radius)
        F = np.cumsum(face_hist)
    else:
        F = 0
    eulers[1:] = V - E + F
    if ec_type == "ECT":
        return eulers
    elif ec_type == "DECT":
        eulers[1:] = (eulers[1:]-eulers[:-1])/(radius[1:]-radius[:-1])
        return eulers
    elif ec_type == "SECT":
        eulers -= np.mean(eulers[i])
        eulers = np.cumsum(eulers)*((radius[-1]-radius[0])/n_filtration)
        return eulers
    else:
        return None

def compute_ec_curve(mesh, directions, n_filtration = 25, ball_radius = 1.0, ec_type = "ECT", first_column_index = False, include_faces = True):
    eulers = np.zeros((directions.shape[0],n_filtration),dtype=float)
    radius = np.linspace(-ball_radius,ball_radius,n_filtration)    
    for i in range(directions.shape[0]):
        eulers[i] = compute_ec_curve_single(mesh,directions[i],radius,n_filtration,ec_type,include_faces)
    return radius, eulers

def compute_ec_curve_parallel(mesh, directions, n_filtration = 25, ball_radius = 1.0, ec_type = "ECT", include_faces = True, n_core = -1):
    radius = np.linspace(-ball_radius,ball_radius,n_filtration)
    parameter = (n_filtration,ball_radius,ec_type,include_faces)
    if n_core == -1:    
        n_core = multiprocessing.cpu_count()
    processed_list = Parallel(n_jobs=n_core)(delayed(compute_ec_curve_single)(mesh,direction,radius,*parameter) for direction in directions)
    processed_list = np.array(processed_list)
    return radius, processed_list

def compute_ec_curve_folder(protA, protB, directions, n_sample = 101, ec_type = "ECT", n_cone = 15, directions_per_cone = 4, cap_radius = 0.1, n_filtration = 25, ball_radius = 1.0, include_faces = True, directory_mesh_A = None, directory_mesh_B = None, sm_radius = 4.0, hemisphere=False, parallel = False, n_core = -1):

    if directory_mesh_A == None or directory_mesh_B == None:
        directory = "%s_%s"
        directory_mesh_A = "%s_%s/mesh/%s_%.1f"%(protA,protB,protA,sm_radius)
        directory_mesh_B = "%s_%s/mesh/%s_%.1f"%(protA,protB,protB,sm_radius) 
        if not os.path.exists(directory):
            os.mkdir(directory)
        if not os.path.exists(directory_mesh_A):
            os.mkdir(directory_mesh_A)
        if not os.path.exists(directory_mesh_B):
            os.mkdir(directory_mesh_B)

        outfiles = []
        for prot, directory_mesh in zip([protA,protB],[directory_mesh_A,directory_mesh_B]):
            ecs = []
            for i_sample in range(n_sample):
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
            if include_faces:
                outfile = '%s/%s_%s_%d_%d_%.1f_%d_%d.txt'%(directory,ec_type,prot,n_cone,directions_per_cone,cap_radius,n_filtration)
            else:
                outfile = '%s/%s_%s_%d_%d_%.1f_%d_%d_nofaces.txt'%(directory,ec_type,prot,n_cone,directions_per_cone,cap_radius,n_filtration)
            
            np.savetxt(outfile,ecs,fmt='%.3f')
            outfiles.append(outfile)

    else:
        directory = '.'
        outfiles = []
        for prot, directory_mesh in zip([protA,protB],[directory_mesh_A,directory_mesh_B]):
            ecs = []
            for filename in os.listdir(directory_mesh):
                if filename.endswith(".msh"):
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
            if include_faces:
                outfile = '%s_%s_%d_%d_%.1f_%d_%d.txt'%(ec_type,prot,n_cone,directions_per_cone,cap_radius,n_filtration)
            else:
                outfile = '%s_%s_%d_%d_%.1f_%d_%d_nofaces.txt'%(ec_type,prot,n_cone,directions_per_cone,cap_radius,n_filtration)

            np.savetxt(outfile,ecs,fmt='%.3f')
            outfiles.append(outfile)

    data_A = np.loadtxt(outfiles[0])
    data_B = np.loadtxt(outfiles[1])

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
    
    normec_file = '%s/%s_%s_%s_%d_%d_%.1f_%d_%d_norm.txt'%(directory,ec_type,protA,protB,n_cone,directions_per_cone,cap_radius,n_filtration)
    notvacuum_file = '%s/notvacuum_%s_%s_%s_%d_%d_%.1f_%d_%d_norm.txt'%(directory,ec_type,protA,protB,n_cone,directions_per_cone,cap_radius,n_filtration)
    np.savetxt(normec_file,data)
    np.savetxt(notvacuum_file,not_vacuum,fmt="%d")
     
    n_A = data_A.shape[0]
    n_B = data_B.shape[0]
    label = np.zeros(n_A+n_B,dtype=int)
    label[:n_A].fill(-1)
    label[n_A:].fill(1)

    return data, label, not_vacuum



