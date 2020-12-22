#!/bin/python3

# calculation libraries
import numpy as np
from scipy.spatial import distance

# Visualization libraries
import matplotlib
from matplotlib import rc
from matplotlib import pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

# distance based filtration
class ComplexFiltration:
    
    def __init__(self):
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
    
    # convert edge list to face list
    def edge_to_face_list(self,edges):
        connections = [[] for i in range(self.n_vertices)]
        for edge in edges:
            connections[edge[0]].append(edge[1])
            connections[edge[1]].append(edge[0])     
#        for i in range(self.n_vertices):
#            connections[i] = np.sort(connections[i])
        faces = []
        for i in range(self.n_vertices):
            for a in connections[i]:
                if a < i:
                    continue
                for b in connections[a]:
                    if b < a:
                        continue
                    elif i in connections[b]:
                        faces.append([a,b,i])
        faces = np.array(faces)
        return faces
    
    # output OFF file for visualization
    def write_mesh_file(self,edges,faces,filename='output.mesh'):
        with open(filename,'w') as f:
            f.write('%d %d %d\n'%(self.vertices.shape[0],edges.shape[0],faces.shape[0]))
            for vertex in self.vertices:
                f.write('%.6f %.6f %.6f\n'%(vertex[0],vertex[1],vertex[2]))
            for edge in edges:
                f.write('%d %d\n'%(edge[0],edge[1]))
            for face in faces:
                f.write('%d  %d %d %d  \n'%(len(face),face[0],face[1],face[2]))
        return
   
    # output OFF file for visualization
    def write_off_file(self,edges,faces,filename='output.off'):
        with open(filename,'w') as f:
            f.write('OFF\n')
            f.write('%d %d %d\n'%(self.vertices.shape[0],faces.shape[0],edges.shape[0]))
            for vertex in self.vertices:
                f.write('%.6f %.6f %.6f\n'%(vertex[0],vertex[1],vertex[2]))
            for face in faces:
                f.write('%d  %d %d %d  \n'%(len(face),face[0],face[1],face[2]))
        return

    def plot_EC_curve(self,radii=np.linspace(0,1,21),outfile='ec_filtration.pdf'):
        ECs = []
        for radius in radii:
            edges, distances = self.get_edge_list(radius=radius)
            faces = self.edge_to_face_list(edges=edges)
            V = self.n_vertices
            E = edges.shape[0]
            F = faces.shape[0]
            EC = V - E + F
            ECs.append(EC)
            print('%.1f %d %d %d %d'%(radius,V,E,F,EC))

        fig, ax = plt.subplots()
        ax.plot(radii,ECs,'ko-')
        ax.set_xlim(radii[0],radii[-1])
        ax.set_ylim(0,1000)
        ax.set_xlabel('radius',fontsize=20)
        ax.set_ylabel('Euler characteristics',fontsize=20)
        ax.tick_params(axis='both',labelsize=20)
        plt.tight_layout()
        plt.savefig(outfile)
        plt.show()
        plt.close()
        return

    def plot_faces_2D(self,edges,faces,radius=1.0,filename='output.pdf'):

        fig, ax = plt.subplots()

        vertices = self.vertices[:,:2]
        ax.scatter(vertices[:,0],vertices[:,1],color='k')
        for edge in edges:
            x = [vertices[edge[0],0],vertices[edge[1],0]]
            y = [vertices[edge[0],1],vertices[edge[1],1]]
            ax.plot(x,y,'k')

        patches = []
        for face in faces:
            polygon = Polygon([vertices[face[0]],vertices[face[1]],vertices[face[2]]],True)
            patches.append(polygon)

        p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
        colors = 100*np.random.rand(len(patches))
        p.set_array(np.array(colors))
        ax.add_collection(p)

        ax.set_xlim(-radius,radius)
        ax.set_ylim(-radius,radius)
        ax.set_xlabel(r'$x$',fontsize=20)
        ax.set_ylabel(r'$y$',fontsize=20)
        ax.tick_params(axis='both',labelsize=20)
        plt.tight_layout()
        plt.savefig(filename)
        plt.show()
        plt.close()
        return
    
    def normalize(self):
        rmax = np.amax(np.linalg.norm(self.vertices,axis=1))
        self.vertices /= rmax
        return

    def calc_radius(self):
        return np.amax(np.linalg.norm(self.vertices,axis=1))
