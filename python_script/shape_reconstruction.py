#!/bin/python3

import numpy as np

#### Shape Reconstruction ####

def summarize_vertices(directions,mesh,rate_vals,length,cone_size,threshold=-1,ball=True,ball_radius=1): 
    if threshold == -1:
        threshold = 1./len(rate_vals)  
    indices = np.arange(len(rate_vals))[rate_vals > threshold]
    selected_vertices = np.array([])
    for i in range(directions.shape[0]):
        vtx_projection = np.dot(mesh.vertices,directions[i])
        if ball:
            binsize = 2*ball_radius/length
            projection_bucket = np.floor((vtx_projection+ball_radius)/binsize).astype(int)
        else:
            binmin = np.amin(vtx_projection)
            binsize = (np.amax(vtx_projection)-binmin)/length
            projection_bucket = np.floor((vtx_projection-binmin)/binsize).astype(int)
        projection_bucket += i*length
        selected_vertices = np.append(selected_vertices,np.intersect1d(projection_bucket,indices))
    final_selected_vertices = np.unique(selected_vertices)
    return final_selected_vertices

def compute_selected_vertices_cones(directions,mesh,rate_vals,length=100,threshold=-1,cone_size=10,ball=True,ball_radius=1.0):
    if threshold == -1:
        threshold = 1./len(rate_vals)
    if directions.shape[0] % cone_size != 0:
        print('Number of Cones not a multiple of directions')
        return 0
    coned_vertices = np.array([])
    for j in range(int(directions.shape[0] / cone_size)):
        cone_dirs = directions[j*cone_size:(j+1)*cone_size]
        cone_rate_vals = rate_vals[j*cone_size*length:(j+1)*cone_size*length]
        coned_vertex = summarize_vertices(directions=cone_dirs,mesh=mesh,rate_vals=cone_rate_vals,length=length,threshold=threshold,cone_size=cone_size,ball=ball,ball_radius=ball_radius,radius=radius)
        print(coned_vertex)
        coned_vertices = np.append(coned_vertices,coned_vertex)
    total_selected_vertices=np.unique(coned_vertices)
    return total_selected_vertices

#### Heatmap Code ####
def reconstruct_vertices_on_shape(directions,mesh,rate_vals,length,cuts,cone_size,ball_radius,ball=True,radius=0): 
    return


