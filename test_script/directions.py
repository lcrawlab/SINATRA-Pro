#!/bin/python3

import numpy as np

#  Compute the directions about a cone using the rodrigues angle formula.
#   @export
#  @param z (vector of length 3) The direction around which we want to rotate.
#  @param r (float): The radius of the cone we want. Controls the cone size.
#  @param j (int): The number of directions in the cone.
#  @return dir (jx3 matrix): of directions in a cone of radius r about z.
def rodrigues(z,N):
    directions = []
    A = np.array([[ 0   ,-z[2], z[1]],
                  [ z[2],    0, z[1]],
                  [-z[1], z[0],    0]])
    thetas = np.linspace(0,2*np.pi,N)
    for theta in thetas:
        R = np.identity(3) + A*np.sin(theta) + np.matmul(A,A)*(1-np.cos(theta))
        direction = np.dot(A,z)
        directions.append(direction)
    return directions

#  Generate Equidistributed points on a sphere
# 
#  @export
#  @description Generate equidistributed points about a sphere using (insert reference)
# 
#  @param desired_number (int): the desired number of equidistributed points on the 2-sphere.
#  @param N (int): the initial number of points that the algorithm will try to generate. If the number of points generated is less than the desired number, the function will increment N
#  @return points (desired_number x 3 matrix): The matrix of equidistributed points in on the sphere.
def generate_equidistributed_points(desired_number, N):
    a = 4*np.pi/N
    d = np.sqrt(a)
    points = [] #np.zeros((n_directions,3),dtype=float)
    M_theta = int(round(np.pi/d))
    d_theta = np.pi/M_theta
    d_phi = a/d_theta
    for i in range(M_theta):
        theta = np.pi * (i + 0.5) / M_theta
        M_phi = int(round(2*np.pi*np.sin(theta)/d_phi))
        for j in range(M_phi):
            phi = 2*np.pi*j/M_phi
            point = [np.sin(theta)*np.cos(phi),np.sin(theta)*np.sin(phi),np.cos(theta)]
            point /= np.linalg.norm(point)
            points.append(point)
    points = np.array(points)
    if points.shape[0] < desired_number:
        return generate_equidistributed_points(desired_number,N+1)
    else:
        return points

# Generate Equidistributed Cones
# 
#  @export
# @description  Generate cones that are equidistributed about the sphere. We first generate equidistributed points, then use the rodrigues angle formula
# to compute cones about those points.
# 
# @param num_directions (int): The number of equidistributed directions we want on the sphere.
# @param cap_radius (float): The radius of the cones we generate (determines the size of each cone).
# @param directions_per_cone (int): The number of directions we want generated within each cone.
# 
# @return directions (num_directions*directions_per_cone x 3 matrix): A matrix of equidistributed cones.
def generate_equidistributed_cones(n_directions, cap_radius, directions_per_cone):
    sphere = generate_equidistributed_points(n_directions, n_directions)
    directions = []
    for i in range(n_directions):
        directions.append(sphere[i]*(cap_radius/np.linalg.norm(sphere[i])))
        if directions_per_cone > 1:
            for direction in rodrigues(sphere[i],directions_per_cone-1):
                direction *= cap_radius/np.linalg.norm(direction)
                directions.append(direction)
    return directions


