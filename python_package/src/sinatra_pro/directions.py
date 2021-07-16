#!/bin/python3

import numpy as np

def rodrigues(z,r,j):
    """
    Compute the directions about a cone using the rodrigues angle formula.

   `z` is the vector of the central axis of the cone.

   `r` is the radius of the cone that controls the size of the cone.

   `j` is the number of directions in the cone.

    """
    z /= np.linalg.norm(z)
    z0 = np.zeros(3)
    if np.any(z == 0):
        z0[z == 0] = 1
    else:
        z0[0] = 2/z[0]
        z0[1] = -1/z[1]
        z0[2] = -1/z[2]
    z0 /= np.linalg.norm(z0)
    z0 *= r
    z0 = z + z0
    z0 /= np.linalg.norm(z0)
    B = np.cross(z,z0)
    C = np.dot(z,z0)*z
    directions = np.zeros((j,3))
    for i in range(j):
        x = 2*np.pi*(i+1)/j
        directions[i,:] = z0*np.cos(x)+B*np.sin(x)+C*(1-np.cos(x))
    return directions

def generate_equidistributed_points(desired_number, N, hemisphere=False):
    """
    Generate Equidistributed points on a sphere / hemi-sphere.
    
    `desired_number` is the desired number of equidistributed points on the 2-sphere.
    
    `N` is the initial number of points that the algorithm will try to generate. If the number of points generated is less than the desired number, the function will increment `N`.
    
    If `hemisphere` is set to true, it generates points over hemisphere instead of over whole sphere.
    """
    if hemisphere:
        a = 2*np.pi/N
        d = np.sqrt(a)
        M_theta = int(round(np.pi*.5/d))
        d_theta = np.pi*.5/M_theta
    else:
        a = 4*np.pi/N
        d = np.sqrt(a)
        M_theta = int(round(np.pi/d))
        d_theta = np.pi/M_theta
    d_phi = a/d_theta
    points = []
    for i in range(M_theta):
        if hemisphere:
            theta = np.pi * .5 * i / M_theta
        else:
            theta = np.pi * (i + 0.5) / M_theta
        M_phi = int(round(2*np.pi*np.sin(theta)/d_phi))
        for j in range(M_phi):
            phi = 2*np.pi*j/M_phi
            point = [np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)]
            point /= np.linalg.norm(point)
            points.append(point)
    points = np.array(points)
    if points.shape[0] < desired_number:
        return generate_equidistributed_points(desired_number,N+1,hemisphere)
    else:
        return points

def generate_equidistributed_cones(n_cone, cap_radius = 0.1, n_direction_per_cone = 1, hemisphere=False):
    """
    Generate equidistributed cones that are equidistributed about the sphere. 
    
    It first generates equidistributed points, then use the rodrigues angle formula to compute cones about those points.
    
    `num_directions` is the number of equidistributed directions we want on the sphere.

    `cap_radius` is the radius of the cones we generate (determines the size of each cone).

    `directions_per_cone` is the number of directions we want generated within each cone.
    
    If `hemisphere` is set to True, it generates points over hemisphere instead of over whole sphere.
    """
    sphere = generate_equidistributed_points(n_cone, n_cone, hemisphere)
    directions = []
    for i in range(n_cone):
        directions.append(sphere[i]/np.linalg.norm(sphere[i]))
        if n_direction_per_cone > 1:
            for direction in rodrigues(sphere[i],cap_radius,n_direction_per_cone-1):
                direction /= np.linalg.norm(direction)
                directions.append(direction)
    return np.array(directions)


