#!/bin/python3

import numpy as np
from itertools import combinations

def smaller_in_front(m,n):
    if m > n:
        return n, m
    else:
        return m, n

class mesh:
    def __init__(self):
        self.vertices = []
        self.edges = []
        self.faces = []
        return
    ''' 
    def intersecting_edges(self):

        ## I connect all possible pairs of points on each face
        ## i.e. for squares, I also connect the diagonals
        ## then count number of connecting lines between two point 
        ## since only the "edge" lines will border another face,
        ## the bordering edge will have 2 "borders"
        ## while the diagonals will have 1 "border" and get filtered out
        borders = np.zeros((self.n_vertex,self.n_vertex),dtype=int)
        for face in self.faces:
            for a in combinations(face,2):
                if a[0] < a[1]:
                    borders[a[0],a[1]] += 1
                else:
                    borders[a[1],a[0]] += 1
        
        ## when the connecting line is shared by two faces, it counts as an edge bordering two faces
        ## assuming the mesh does not contain face that, e.g. hangs on the diagonals of a square face
        edges = np.transpose((borders == 2).nonzero())
        n_edges = edges.shape[0]
        
        ## check if all vertices are connected (not isolated)
        print((np.any(borders,axis=1)==0).nonzero())

        ## check if number of edges on the face match the number of vertices (e.g. 3 for a triangle, 4 for a sqaure)
        ## if the mesh is continuous, there will be no mismatch
        mismatch = False
        for face in self.faces:
            n_edgeConnected = 0
            for a in combinations(face,2):
                if (a[0] < a[1] and borders[a[0]-1,a[1]-1] == 2) or (a[0] > a[1] and borders[a[1]-1,a[0]-1] == 2):
                    n_edgeConnected += 1
            if n_edgeConnected != len(face):
                n_edges += len(face) - n_edgeConnected
                mismatch = True
        if mismatch:
            print("There are mismatch between number of edges and number of edges that connect two faces. The mesh might be unclosed.")
        else:
            print("There are no mismatch")
        
        ##
        ## Edges are defined as connection between any two points shared by two faces
        ##
        borders = []
        edges = []
        for face in self.faces:
            for a in combinations(face,2):
                m, n = smaller_in_front(a[0],a[1])
                if [m,n] in borders:
                    edges.append([m,n])
                    borders.remove([m,n])
                else:
                    borders.append([m,n])
        
        ## 
        ## check if number of edges on the face match the number of vertices (e.g. 3 for a triangle, 4 for a sqaure)
        ## if the mesh is continuous, there will be no mismatch
        ##
        mismatch = False
        missed_edges = 0
        for face in self.faces:
            n_edgeConnected = 0
            for a in combinations(face,2):
                m, n = smaller_in_front(a[0],a[1])
                if [m,n] in edges:
                    n_edgeConnected += 1
            if n_edgeConnected != len(face):
                missed_edges += len(face) - n_edgeConnected
                mismatch = True
        edges = np.array(edges)
        n_edges = edges.shape[0]
        if mismatch:
            print("There are mismatch between number of edges and number of edges that connect two faces. The mesh might be unclosed.")
        else:
            print("There are no mismatch. The mesh is closed.")
        return edges
    '''
    def unique_edges(self):
        edges = []
        for face in self.faces:
            for i,j in [[0,1],[1,2],[2,3],[3,0]]:
                m, n = smaller_in_front(face[i],face[j])
                edges.append([m,n])
        self.edges = np.unique(edges,axis=0)
        return

    def read_off_file(self,filename,read_edges_from_file=False,edges_filename='edges.txt'):
        vertices = []
        faces = []
        nline = 0
        n_V = 0
        n_F = 0
        n_E = 0
        with open(filename,'r') as f:
            for line in f:
                if line[0] == '#' or line.strip() == "OFF":
                    continue
                p = list(filter(lambda x: x != '',line.strip().split(' ')))
                if nline == 0:
                    n_V = int(p[0])
                    n_F = int(p[1])
                    n_E = int(p[2])
                elif nline < n_V + 1:
                    vertex = np.array([float(a) for a in p])
                    vertices.append(vertex)
                else:
                    #if int(p[0]) > 3:
                    #    nontriangular = True
                    face = [int(p[i+1]) for i in range(int(p[0]))]
                    faces.append(face)
                nline += 1
        self.n_vertex = n_V
        self.n_face = n_F
        self.vertices = np.array(vertices)
        self.faces = np.array(faces)
        if read_edges_from_file:
            self.edges = np.loadtxt(edges_filename,dtype=int,delimiter=' ')
        else:
            self.unique_edges()
        self.n_edge = self.edges.shape[0]
        return
   
    def read_mesh_file(self,filename='output.mesh'):
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
                    self.edges = np.zeros((self.n_edge,2),dtype=int)
                    self.faces = []
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
                        face = []
                        for i in range(int(p[0])):
                            face.append(int(p[i+1]))
                        self.faces.append(face)
                n_line += 1
        return
    
    def write_off_file(self,filename='output.off'):
        with open(filename,'w') as f:
            f.write('OFF\n')
            f.write('%d %d 0\n'%(self.vertices.shape[0],self.faces.shape[0]))
            for vertex in self.vertices:
                f.write('%.6f %.6f %.6f\n'%(vertex[0],vertex[1],vertex[2]))
            for face in self.faces:
                f.write('%d  '%len(face))
                for a in face:
                    f.write('%d '%a)
                f.write(' \n')
        return
    
    def write_mesh_file(self,edges,faces,filename='output.mesh'):
        with open(filename,'w') as f:
            f.write('%d %d %d\n'%(self.vertices.shape[0],edges.shape[0],faces.shape[0]))
            for vertex in self.vertices:
                f.write('%.6f %.6f %.6f\n'%(vertex[0],vertex[1],vertex[2]))
            for edge in edges:
                f.write('%d %d\n'%(edge[0],edge[1]))
            for face in faces:
                f.write('%d  '%len(face))
                for a in face:
                    f.write('%d '%a)
                f.write(' \n')
        return
   
    def centering(self):
        center = np.mean(self.vertices,axis=0)
        self.vertices = self.vertices - center
        return
    
    def write_off_file_heat(self,heat,filename='output.off'):
        with open(filename,'w') as f:
            f.write('OFF\n')
            f.write('%d %d 0\n'%(self.vertices.shape[0],self.faces.shape[0]))
            for vertex in self.vertices:
                f.write('%.6f %.6f %.6f\n'%(vertex[0],vertex[1],vertex[2]))
            for i, face in enumerate(self.faces):
                f.write('%d  '%len(face))
                for a in face:
                    f.write('%d '%a)
                f.write(' %d %d %d \n'%(heat[i,0],heat[i,1],heat[i,2]))
        return

