#!/bin/python3

import numpy as np

def compute_vertex_function(mesh,direction):
    vf = np.dot(mesh.vertices,direction)
    return vf

#' Computes the EC curve on a ball
#'
#' @export
#' @description Computes the EC curve of a mesh on a bounding ball of specified radius $radius$. The filtration steps for which the
#' (S/D) EC curve are computed is done relative to the bounding ball. For filtrations steps that are computed relative to the shape, use
#' \code{compute_discrete_ec_curve}. For comparisons with multiple shapes, this is the recommended function.
#'
#' @param complex (list) : The list containing metadata about the mesh; use \code{process_off_filev3} to obtain this list from an off file, or
#' \code{convert_complex} to convert existing mesh3d objects to this list.
#' @param vertex_function (matrix): A matrix containing the projections of each vertex along a given direction. Computed as the dot product of the vertex onto the
#' direction on the unit sphere.
#' @param curve_length (int): The number of sub-level sets for which to compute the EC curve on in a given direction.
#' @param first_column_index (boolean): Specifying the vertex index is included in the vertex function matrix.
#' @param ball_radius (float): The radius of the bounding ball.
#'
#' @return ec_curve (nx2 matrix): A matrix containing the EC curve of the given simplicial complex with the first index as the projections for which the EC was computed on.

def compute_ec_curve(mesh, vertex_function, curve_length = 20, ball_radius = 2.0, ec_type = "EC", standardized = False,first_column_index = False):
    if standardized:
        radius = np.linspace(-ball_radius,ball_radius,curve_length+1)
    else:
        radius = np.linspace(np.amin(vertex_function),np.amax(vertex_function),curve_length+1)
    vertex_hist, bin_edges = np.histogram(vertex_function,bins=radius)
    V = np.cumsum(vertex_hist)
    # filtrating edges
    edge_function = np.amax(vertex_function[mesh.edges],axis=1)
    edge_hist, bin_edges = np.histogram(edge_function,bins=radius)
    E = np.cumsum(edge_hist)
    # filtrating faces
    face_function = np.array([np.amax(vertex_function[face]) for face in mesh.faces])
    face_hist, bin_edges = np.histogram(face_function,bins=radius)
    F = np.cumsum(face_hist)
    euler = V-E+F
    if ec_type == "EC":
        return radius, euler
    elif ec_type == "DECT":
        return radius, np.append(0,differentiate(x,euler))
    elif ec_type == "SECT":
        return radius, np.appedn(0,integrate(x,euler))
    else:
        return

def differentiate(x,y):
    return (y[1:]-y[:-1])/(x[1:]-x[:-1])

def integrate(x,y):
    trapz = (y[1:]+y[:-1])*(x[1:]+x[:-1])*.5
    return np.cumsum(trapz)/y.size
