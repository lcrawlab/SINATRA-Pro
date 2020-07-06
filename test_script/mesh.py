#!/bin/python3

import numpy as np
from itertools import combinations


class mesh:
    def __init(self):
        return

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
                    borders[a[0]-1,a[1]-1] += 1
                else:
                    borders[a[1]-1,a[0]-1] += 1
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
        self.faces = faces #np.array(faces)
        self.edges = self.intersecting_edges()
        return
