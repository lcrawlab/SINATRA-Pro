#!/bin/python3

import numpy as np
from scipy.integrate import cumtrapz 

#' Computes the EC curve on a ball
#'
#' @export
#' @description Computes the EC curve of a mesh on a bounding ball of specified radius $radius$. The filtration steps for which the
#' (S/D) EC curve are computed is done relative to the bounding ball. For filtrations steps that are computed relative to the shape, use
#' \code{compute_discrete_ec_curve}. For comparisons with multiple shapes, this is the recommended function.
#'
#' @param complex (list) : The list containing metadata about the mesh; use \code{process_off_filev3} to obtain this list from an off file, or
#' \code{convert_complex} to convert existing mesh3d objects to this list.
#' @param direction (array)
#####' @param vertex_function (matrix): A matrix containing the projections of each vertex along a given direction. Computed as the dot product of the vertex onto the
#' direction on the unit sphere.
#' @param curve_length (int): The number of sub-level sets for which to compute the EC curve on in a given direction.
#' @param ball_radius (float): The radius of the bounding ball.
#'
#' @return ec_curve (nx2 matrix): A matrix containing the EC curve of the given simplicial complex with the first index as the projections for which the EC was computed on.


def compute_ec_curve(mesh, directions, n_filtration = 20, ball_radius = 2.0, ec_type = "EC", standardized = False,first_column_index = False):

    n_direction = directions.shape[0]
    eulers = np.zeros((n_direction,n_filtration+1),dtype=float)
     
    for i in range(n_direction): 
        vertex_function = np.dot(mesh.vertices,directions[i]) 
        if standardized:
            radius = np.linspace(-ball_radius,ball_radius,n_filtration+1)
        else:
            radius = np.linspace(np.amin(vertex_function),np.amax(vertex_function),n_filtration+1)
        # filtrating vertices
        vertex_hist, bin_edges = np.histogram(vertex_function,bins=radius)
        V = np.cumsum(vertex_hist)
        # filtrating edges
        edge_function = np.amax(vertex_function[mesh.edges],axis=1)
        edge_hist, bin_edges = np.histogram(edge_function,bins=radius)
        E = np.cumsum(edge_hist)
        # filtrating faces
        face_function = [np.amax(vertex_function[face]) for face in mesh.faces]
        face_hist, bin_edges = np.histogram(face_function,bins=radius)
        F = np.cumsum(face_hist)
        eulers[i] = np.append(0,V-E+F)
    
    if ec_type == "EC":
        return radius, eulers
    elif ec_type == "DECT":
        dect = np.zeros((n_direction,n_filtration+1),dtype=float)
        for i in range(n_direction):
            dect[i,1:] = (eulers[i,1:]-eulers[i,:-1])/(radius[1:]-radius[:-1])
        return radius, dect
    elif ec_type == "SECT":
        sect = np.zeros((n_direction,n_filtration),dtype=float)
        for i in range(n_direction):
            sect[i] = cumtrapz(euler,radius,initial=0)
        return radius, sect
    else:
        return

