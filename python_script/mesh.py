#!/bin/python3

import numpy as np

class mesh:
    def __init(self):
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
    
    def centering(self):
        center = np.mean(self.vertices,axis=0)
        self.vertices = self.vertices - center
        return
