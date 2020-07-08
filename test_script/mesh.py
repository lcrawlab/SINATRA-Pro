#!/bin/python3

import numpy as np
from itertools import combinations

import timeit

def smaller_in_front(m,n):
    if m > n:
        return n, m
    else:
        return m, n

class mesh:
    def __init(self):
        return

    ## new version to save memory for large mesh
    def intersecting_edges(self):
        ## I connect all possible pairs of points on each face
        ## i.e. for squares, I also connect the diagonals
        ## then count number of connecting lines between two point 
        ## since only the "edge" lines will border another face,
        ## the bordering edge will have 2 "borders"
        ## while the diagonals will have 1 "border"
        borders = []
        edges = []
        for face in self.faces:
            for a in combinations(face,2):
                m, n = smaller_in_front(a[0],a[1])
                #print(borders,[m,n],[m,n] in borders)
                if [m,n] in borders:
                    edges.append([m,n])
                    borders.remove([m,n])
                else:
                    borders.append([m,n])

        ## when the connecting line is shared by two faces, it counts as an edge bordering two faces
        ## assuming the mesh does not contain face that, e.g. hangs on the diagonals of a square face
        #edges = np.transpose((borders == 2).nonzero())

        edges = np.array(edges)
        print(borders)
        print(edges.shape)
        exit()
        n_edges = edges.shape[0]
        
        ## check if all vertices are connected (not isolated)
        print(borders)

        ## check if number of edges on the face match the number of vertices (e.g. 3 for a triangle, 4 for a sqaure)
        ## if the mesh is continuous, there will be no mismatch
        mismatch = False
        for face in self.faces:
            n_edgeConnected = 0
            for a in combinations(face,2):
                m, n = smaller_in_front(a[0],a[1])
                if np.any(edges == [m,n]):
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

    def read_off_file(self,filename):
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
                    face = [int(p[i+1]) for i in range(int(p[0]))]
                    faces.append(face)
                nline += 1
        self.n_vertex = n_V
        self.n_face = n_F
        self.n_edge = n_E
        self.vertices = np.array(vertices)
        self.faces = np.array(faces)
        self.edges = self.intersecting_edges()
        return

    def read_obj_file(self,filename,nontriangular=False):
        vertices = []
        faces = []
        nline = 0
        with open(filename,'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                p = list(filter(lambda x: x != '',line.strip().split(' ')))
                if p[0] == 'v':
                    vertex = [float(a) for a in p[1:]]
                    vertices.append(vertex)
                elif p[0] == 'f': 
                    face = [int(a.split('/')[0]) for a in p[1:]]
                    faces.append(face)
                nline += 1
        self.vertices = np.array(vertices)
        self.faces = np.array(faces)
        if nontriangular == True:
            self.edges = self.intersecting_edges()
        else:
            edges = []
            for face in self.faces:
                for a in combinations(face,2):
                    m, n = smaller_in_front(a[0],a[1])
                    edges.append([m,n])
            self.edges = np.unique(edges,axis=0)
        self.n_vertex = self.vertices.shape[0]
        self.n_face = self.faces.shape[0]
        self.n_edge = self.edges.shape[0]
        return


