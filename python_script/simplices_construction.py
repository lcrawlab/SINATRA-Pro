#!/bin/python3

# calculation libraries
import numpy as np
from scipy.spatial import distance

# Visualization libraries
#import matplotlib
#from matplotlib import rc
#from matplotlib import pyplot as plt
#from matplotlib.patches import Polygon
#from matplotlib.collections import PatchCollection

from MDAnalysis.lib.nsgrid import FastNS, NSResults

import multiprocessing
from joblib import Parallel, delayed

# distance based filtration
class ComplexFiltration:
    
    def __init__(self):
        self.vertices = None
        self.n_vertices = None
        self.faces = None
        return
    
    # read positions of vertices from file
    def import_vertices_from_file(self, filename):
        self.vertices = np.loadtxt(filename,usecols=(0,1,2))
        self.n_vertices = self.vertices.shape[0]
        return
    
    # generate N random vertices on [0,1],[0,1]
    def generate_random_vertices(self, n_vertices):
        np.random.seed(0)
        self.n_vertices = n_vertices
        vertices = []
        for i in range(n_vertices):
            vertice = []
            for k in range(3):
                vertice.append(np.random.random())
            vertices.append(vertice)
        vertices = np.array(vertices)
        self.vertices = vertices
        return

    # calculate distance matrix
    def calc_distance_matrix(self):
        self.n_vertices = self.vertices.shape[0]

        # calculate self distance matrix among vertices
        self.distance_matrix = distance.pdist(self.vertices)
        
        # generate pair list that matches the format of the distance matrix from scipy
        # i.e. [[0,1],[0,2] ....,[0,N],[1,2],...,[N-2,N-1]]
        pairs = []
        for i in range(self.n_vertices-1):
            for j in range(i+1,self.n_vertices):
                pairs.append([i,j])
        self.pairs = np.array(pairs)
        
        # sort both pair list and distance matrix by distance
        sorted_filter = np.argsort(self.distance_matrix)
        self.distance_matrix = self.distance_matrix[sorted_filter]
        self.pairs = self.pairs[sorted_filter]
        
        return
    
    # output list of edges within distance cutoff from sorted distance matrix  
    def get_edge_list(self, radius):
        cutoff = np.argmax(self.distance_matrix > radius)
        return self.pairs[:cutoff], self.distance_matrix[:cutoff]
    
    def neighbor_search(self,cutoff=4.0):
        max_gridsize = 8000
        passed = False
        while not passed:
            try:
                nsr = FastNS(cutoff=cutoff,coords=self.vertices,box=self.box,pbc=False,max_gridsize=max_gridsize)
                passed = True
            except MemoryError:
                max_gridsize = int(max_gridsize/2)
            except:
                raise
            if passed:
                break
        result = nsr.self_search()
        self.edges = result.get_pairs()[::2]
        return 
    
    # convert edge list to face list
    def edge_to_face_list(self):
        self.n_vertices = self.vertices.shape[0]
        self.connections = [set() for i in range(self.n_vertices)]
        for edge in self.edges:
            self.connections[edge[0]].add(edge[1])
            self.connections[edge[1]].add(edge[0])
        self.faces = []
        self.checked = set()
        for u in range(self.n_vertices):
            self.checked.add(u)
            for v in self.connections[u] - self.checked:
                for s in self.connections[u] & self.connections[v]:
                    if s > v and v > u:
                        self.faces.append([u,v,s])
        self.faces = np.array(self.faces)
        return
    
    # output OFF file for visualization
    def write_mesh_file(self,filename='output.mesh'):
        with open(filename,'w') as f:
            f.write('%d %d %d\n'%(self.vertices.shape[0],self.edges.shape[0],self.faces.shape[0]))
            for vertex in self.vertices:
                f.write('%.6f %.6f %.6f\n'%(vertex[0],vertex[1],vertex[2]))
            for edge in self.edges:
                f.write('%d %d\n'%(edge[0],edge[1]))
            for face in self.faces:
                f.write('%d  %d %d %d  \n'%(len(face),face[0],face[1],face[2]))
        return
   
    # output OFF file for visualization
    def write_off_file(self,filename='output.off'):
        with open(filename,'w') as f:
            f.write('OFF\n')
            f.write('%d %d %d\n'%(self.vertices.shape[0],self.faces.shape[0],self.edges.shape[0]))
            for vertex in self.vertices:
                f.write('%.6f %.6f %.6f\n'%(vertex[0],vertex[1],vertex[2]))
            for face in self.faces:
                f.write('%d  %d %d %d  \n'%(len(face),face[0],face[1],face[2]))
        return

    def normalize(self):
        rmax = np.amax(np.linalg.norm(self.vertices,axis=1))
        self.vertices /= rmax
        return

    def calc_radius(self):
        return np.amax(np.linalg.norm(self.vertices,axis=1))
