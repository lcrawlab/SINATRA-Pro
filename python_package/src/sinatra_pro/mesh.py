#!/bin/python3

import numpy as np
from scipy.spatial import distance
import MDAnalysis
from MDAnalysis.lib.nsgrid import FastNS, NSResults
import multiprocessing
from joblib import Parallel, delayed

class mesh:
    
    def __init__(self):
        self.vertices = []
        """List of coordinates of vertices"""
        self.n_vertices = 0
        """Number of vertices"""
        self.edges = []
        """List of indices of connected edges"""
        self.n_edge = 0
        """Number of edges"""
        self.faces = []
        """List of indices of constructed faces"""
        self.n_face = 0
        """Number of faces"""
        return
    
    def calc_radius(self):
        """Calculate distance of the vertex furthest away from origin """
        return np.amax(np.linalg.norm(self.vertices,axis=1))
    
    def normalize(self):
        """Normalize vertices to the unit sphere """
        rmax = self.calc_radius()
        self.vertices /= rmax
        return
   
    def centering(self):
        """Center protein to origin by center of geometry"""
        center = np.mean(self.vertices,axis=0) # center of geometry
        self.vertices = self.vertices - center
        return 
    
    def import_vertices_from_file(self, filename):
        """Read positions of vertices from file"""
        self.vertices = np.loadtxt(filename,usecols=(0,1,2))
        self.n_vertices = self.vertices.shape[0]
        return

    def calc_distance_matrix(self):
        """Calculate distance matrix"""
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
    
    def get_edge_list(self, radius):
        """Output list of edges within distance cutoff from sorted distance matrix"""
        pair_list = np.argmax(self.distance_matrix > radius)
        return self.pairs[:pair_list], self.distance_matrix[:pair_list]
    
    def neighbor_search_old(self,coords,box,cutoff):
        """
        Neighbor search using Neighbor Grid Search (FastNS) algorithm from MDAnalysis (==0.20.1). 
        
        It generates a list of indices of pairs of vertices that are within `cutoff` apart, 
        then save the list as list of edges for the mesh.
        """
        nsr = FastNS(cutoff=cutoff,coords=coords,box=box,pbc=False)
        result = nsr.self_search()
        self.edges = result.get_pairs()[::2]
        return 
    
    def neighbor_search_new(self,coords,box,cutoff):
        """
        Neighbor search using Neighbor Grid Search (FastNS) algorithm from MDAnalysis (==1.1.1). 
        
        It generates a list of indices of pairs of vertices that are within `cutoff` apart, 
        then save the list as list of edges for the mesh.
        """
        nsr = FastNS(cutoff=cutoff,coords=coords,box=box,pbc=False)
        result = nsr.self_search()
        self.edges = result.get_pairs()
        return 
   
    def edge_to_face_list(self):
        """
        Convert edge list to face list

        Iterate through list of connected edges to look for any 3 edges that enclose a triangle. 
        Such enclosed triangles are constructed as faces.
        """
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
     
    def read_mesh_file(self,filename='output.mesh'):
        """Read topology from .msh file"""
        n_line = 0
        i_v = 0
        i_e = 0
        i_f = 0
        with open(filename,'r') as f:
            n_line = 0
            for line in f:
                p = list(filter(lambda x: x != '',line.strip().split(' ')))
                if n_line == 0:
                    self.n_vertex = int(p[0])
                    self.n_edge = int(p[1])
                    self.n_face = int(p[2])
                    self.vertices = np.zeros((self.n_vertex,3),dtype=float)
                    self.faces = np.zeros((self.n_face,3),dtype=int)
                    self.edges = np.zeros((self.n_edge,2),dtype=int)
                elif n_line < self.n_vertex + 1:
                    if i_v < self.n_vertex:
                        for i in range(3):
                            self.vertices[i_v][i] = float(p[i])
                        i_v += 1
                elif n_line < self.n_vertex + self.n_edge + 1:
                    if i_e < self.n_edge:
                        for i in range(2):
                            self.edges[i_e][i] = int(p[i])
                        i_e += 1
                else:
                    if i_f < self.n_face:
                        for i in range(int(p[0])):
                            self.faces[i_f][i] = int(p[i+1])
                        i_f += 1
                n_line += 1
        return 
   
    def write_mesh_file(self,filename='output.mesh'):
        """Write topology into .msh file"""
        with open(filename,'w') as f:
            f.write('%d %d %d\n'%(self.vertices.shape[0],self.edges.shape[0],self.faces.shape[0]))
            for vertex in self.vertices:
                f.write('%.6f %.6f %.6f\n'%(vertex[0],vertex[1],vertex[2]))
            for edge in self.edges:
                f.write('%d %d\n'%(edge[0],edge[1]))
            for face in self.faces:
                f.write('%d  %d %d %d  \n'%(len(face),face[0],face[1],face[2]))
        return
   
    def write_off_file(self,filename='output.off'):
        """Write topology into .off file for visualization"""
        with open(filename,'w') as f:
            f.write('OFF\n')
            f.write('%d %d %d\n'%(self.vertices.shape[0],self.faces.shape[0],self.edges.shape[0]))
            for vertex in self.vertices:
                f.write('%.6f %.6f %.6f\n'%(vertex[0],vertex[1],vertex[2]))
            for face in self.faces:
                f.write('%d  %d %d %d  \n'%(len(face),face[0],face[1],face[2]))
        return
    
    def convert_vertices_to_mesh(self,sm_radius=2.0,msh_file='mesh.msh',rmax=1.0):
        """
        Convert set of vertices to a simplicial complex by connecting edges and faces
        
        `sm_radius` is the radius cutoff for constructing simplicial complices. 
        Pairs of vertices closer than `sm_radius` Angstrom apart are connected to form edges. 

        `msh_file` is the filename for output .msh files.

        `rmax` is the radius of the largest mesh used to normalize all meshes to the same unit sphere. 
        """
        # position vertices for grid search algorithm
        temp = self.vertices.copy()
        lmax = np.amax(self.vertices,axis=0)
        lmin = np.amin(self.vertices,axis=0)
        box = np.append((lmax-lmin)*1.2,[90.0,90.0,90.0])
        temp -= lmin

        mda_version = [int(a) for a in MDAnalysis.__version__.split('.')]
        if mda_version[0] == 0 and mda_version[1] >= 19:
            self.neighbor_search_old(cutoff=sm_radius,coords=temp,box=box) # neighbor grid search to identify vertex pairs within r < cutoff apart
        elif mda_version[0] == 1 and (mda_version[2] >= 2 or mda_version[1] >= 1):
            self.neighbor_search_new(cutoff=sm_radius,coords=temp,box=box) # neighbor grid search to identify vertex pairs within r < cutoff apart  
        else:
            self.calc_distance_matrix()
            self.edges, distances = self.get_edge_list(radius=sm_radius)
        self.edge_to_face_list() # generate faces enclosed by any 3 edges    
        self.vertices /= rmax # normalized generated meshes to the specified unit sphere
        self.write_mesh_file(filename=msh_file)
        return
    
    def generate_random_vertices(self, n_vertices):
        """Generate N random vertices on [0,1],[0,1] for testing purpose"""
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
