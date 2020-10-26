#!/bin/python3

import numpy as np
from itertools import combinations

def smaller_in_front(m,n):
    if m > n:
        return n, m
    else:
        return m, n

class mesh:
    def __init(self):
        return

    ## new version to save memory for large mesh data
    
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
                print(face,n_edgeConnected)
                mismatch = True
        print(n_edges)
        if mismatch:
            print("There are mismatch between number of edges and number of edges that connect two faces. The mesh might be unclosed.")
        else:
            print("There are no mismatch")
        return edges

    '''
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
            for a in combinations(face,2):
                m, n = smaller_in_front(a[0],a[1])
                edges.append([m,n])
        edges = np.unique(edges,axis=0)
        return edges

    def read_off_file(self,filename,read_edges_from_file=False,edges_filename='edges.txt'):
        vertices = []
        faces = []
        nline = 0
        n_V = 0
        n_F = 0
        n_E = 0
        nontriangular = False
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
                    if int(p[0]) > 3:
                        nontriangular = True
                    face = [int(p[i+1]) for i in range(int(p[0]))]
                    faces.append(face)
                nline += 1
        self.n_vertex = n_V
        self.n_face = n_F
        self.vertices = np.array(vertices)
        self.faces = faces
        if read_edges_from_file:
            self.edges = np.loadtxt(edges_filename,dtype=int,delimiter=' ')
        else:
            if nontriangular:
                self.edges = self.intersecting_edges()
            else:
                self.edges = self.unique_edges()
        self.n_edge = self.edges.shape[0]
        return

    def read_obj_file(self,filename,read_edges_from_file=False,edges_filename='edges.txt'):
        vertices = []
        faces = []
        nline = 0
        nontriangular = False 
        with open(filename,'r') as f:
            for line in f:
                if line[0] == '#' or len(line.strip()) == 0:
                    continue
                p = list(filter(lambda x: x != '',line.strip().split(' ')))
                if p[0] == 'v':
                    vertex = [float(a) for a in p[1:]]
                    vertices.append(vertex)
                elif p[0] == 'f': 
                    face = [int(a.split('/')[0])-1 for a in p[1:]]
                    if len(face) > 3:
                        nontriangular = True
                    faces.append(face)
                nline += 1
        self.vertices = np.array(vertices)
        self.faces = faces
        if read_edges_from_file:
            self.edges = np.loadtxt(edges_filename,dtype=int,delimiter=' ')
        else:
            if nontriangular:
                self.edges = self.intersecting_edges()
            else:
                self.edges = self.unique_edges()
        self.n_vertex = self.vertices.shape[0]
        self.n_face = len(faces)
        self.n_edge = self.edges.shape[0]
        return

    def save_edges_file(self,filename):
        np.savetxt(filename,self.edges,fmt=['%d','%d'])
        return
    
    def centering(self):
        center = np.mean(self.vertices,axis=0)
        self.vertices = self.vertices - center
        return
    
    def normalize(self):
        rmax = np.amax(np.linalg.norm(self.vertices,axis=1))
        self.vertices /= rmax
        return

    def write_off_file(self,filename):
        with open(filename,'w') as f:
            f.write('OFF\n')
            f.write('%d %d %d\n'%(self.n_vertex,self.n_face,self.n_edge))
            for vertex in self.vertices:
                f.write('%.5f %5f %5f\n'%(vertex[0],vertex[1],vertex[2]))
            for face in self.faces:
                f.write('%d %d %d %d\n'%(len(face),face[0],face[1],face[2]))
        return

    def write_off_file_heat(self,filename,vertices_heat):
        with open(filename,'w') as f:
            f.write('OFF\n')
            f.write('%d %d %d\n'%(self.n_vertex,self.n_face,self.n_edge))
            for vertex in self.vertices:
                f.write('%.5f %5f %5f\n'%(vertex[0],vertex[1],vertex[2]))
            for face in self.faces:
                f.write('%d %d %d %d  %d %d %d\n'%(len(face),face[0],face[1],face[2],heat[face[0],0],heat[face[0],1],heat[face[0],2]))
        return

meshProtein = mesh()
meshProtein.read_off_file('off_noh/WT/WT_4_frame10000.off')

import csv

heat = np.zeros((meshProtein.n_vertex,3),dtype=int)
with open("vert_heat_WT_R164G_noh_unit_15_01_4_25_varvs_pip_WT_10000.txt","r") as f:
    p = csv.reader(f)
    i = 0
    for row in p:
        color = np.array([int(row[0][i:i+2],16) for i in (1, 3, 5)])
        heat[i,:] = color
        i += 1
        
meshProtein.write_off_file_heat('WT_R164G_varbvs.off',heat)

    
