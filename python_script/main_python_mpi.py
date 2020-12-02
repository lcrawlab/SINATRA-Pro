#!/bin/python3

from mesh import *
from euler import *
from directions import *

import numpy as np
import os, sys

from mpi4py import MPI

mutant = "R164S"
ec_type = "DECT"

n_file = 101
n_cone = 15
n_filtration = 100
ball_radius = 1.0

comm = MPI.COMM_WORLD
mpi_size = comm.Get_size()
mpi_rank = comm.Get_rank()
print('Running on rank %d of %d...'%(mpi_rank,mpi_size))

radius = None
meshProtein = None
directions = None
ecs = None
n_direction = None
directions_block = None
ecs_block = None

abort = False
if mpi_rank == 0:
    directions = generate_equidistributed_cones(n_directions=n_cone, directions_per_cone=4, cap_radius = 0.1, hemisphere=False)
    np.savetxt('directions_cone_%d.txt'%n_cone,directions)
    directions = np.loadtxt('directions_cone_%d.txt'%n_cone)
    n_direction = directions.shape[0]
    if n_direction % mpi_size != 0:
        print('# of directions not multiple of mpi cores!')
        abort = True
if abort:
    exit()
n_direction = comm.bcast(n_direction,root=0)

if mpi_rank == 0:
    ecss = np.zeros((n_file,n_direction*n_filtration),dtype=float)    

for i in range(n_file):
    if mpi_rank == 0:
        frame = i*100
        print('Calculating %s for frame %d'%(ec_type,frame))
        radius = np.linspace(-ball_radius,ball_radius,n_filtration)
        meshProtein = mesh()
        meshProtein.read_mesh_file(filename='mesh_noh/%s_65_230_4/%s_4_frame%d.msh'%(mutant,mutant,frame)) 
        ecs = np.zeros((n_direction,n_filtration),dtype=float)

    radius = comm.bcast(radius,root=0)
    meshProtein = comm.bcast(meshProtein,root=0)
    
    for j in range(int(n_direction/mpi_size)):
        if mpi_rank == 0:
            directions_block = directions[j*mpi_size:(j+1)*mpi_size,:]
            ecs_block = ecs[j*mpi_size:(j+1)*mpi_size,:]
        directions_block = comm.scatter(directions_block,root=0)
        ecs_block = comm.scatter(ecs_block,root=0)
        ecs_block = compute_ec_curve_single(meshProtein,directions_block,radius,n_filtration=n_filtration,ball_radius=ball_radius,ec_type=ec_type,include_faces=True)
        ecs_block = comm.gather(ecs_block,root=0)
        if mpi_rank == 0:
            ecs[j*mpi_size:(j+1)*mpi_size,:] = ecs_block
 
    if mpi_rank == 0:
        ecss[i] = ecs.flatten()

if mpi_rank == 0:
    np.savetxt('%s_%s_65_230_4_%d_01_4_%d.txt'%(ec_type,mutant,n_cone,n_filtration),ecss,fmt='%.3f')
