#!/bin/python3

import numpy as np
import multiprocessing
from joblib import Parallel, delayed

def compute_ec_curve_single(mesh, direction, radius, n_filtration = 25, ball_radius = 1.0, ec_type = "ECT", include_faces = True):
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
        eulers = (eulers[1:]-eulers[:-1])/(radius[1:]-radius[:-1])
        return eulers
    elif ec_type == "SECT":
        eulers -= np.mean(eulers[i])
        eulers = np.cumsum(eulers)*(radius[1]-radius[0])
        return eulers
    else:
        return None

def compute_ec_curve_parallel(mesh, directions, n_filtration = 25, ball_radius = 1.0, ec_type = "ECT", include_faces = True, n_core = -1):
    radius = np.linspace(-ball_radius,ball_radius,n_filtration)
    parameter = (n_filtration,ball_radius,ec_type,include_faces)
    if n_core == -1:    
        n_core = multiprocessing.cpu_count()
    processed_list = Parallel(n_jobs=n_core)(delayed(compute_ec_curve_single)(mesh,direction,radius,*parameter) for direction in directions)
    processed_list = np.array(processed_list)
    return radius, processed_list


