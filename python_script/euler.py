#!/bin/python3

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
