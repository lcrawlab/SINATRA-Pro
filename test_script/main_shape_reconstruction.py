#!/bin/python3

import numpy as np
from mesh import *
from shape_reconstruction import *

directions = np.loadtxt('directions_cone.txt')
meshProtein = mesh()
meshProtein.read_obj_file(filename='/users/wtang8/scratch/TDA/mesh/rssb/WT/WT_chimera_1.obj')
rate_vals = np.loadtxt('rates.txt',usecols=1)
summarize_vertices(directions=directions,mesh=meshProtein,rate_vals=rate_vals,length=100,cone_size=10,ball=True,ball_radius=50)



