#!/bin/python3

from distance_based_filtration import *
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import numpy as np, os, sys

mutants = ["WT","R164S","R164G"]
nsample = 100
'''
selection = 'resid 65:230'

refmutant = "WT"
refu = mda.Universe('%s/md_0_1.gro'%refmutant,'%s/md_0_1_noPBC.xtc'%refmutant)
refu.trajectory[-1]
refuCA = refu.select_atoms('name CA and %s'%selection)
com_refuCA = refuCA.center_of_mass()
refu0 = refuCA.positions - com_refuCA

for mutant in mutants:

    u = mda.Universe('%s/md_0_1.gro'%mutant,'%s/md_0_1_noPBC.xtc'%mutant)

    u.trajectory[-1]
    refCA = u.select_atoms('name CA and %s'%selection)
    com_refCA = refCA.center_of_mass()
    ref0 = refCA.positions - com_refCA

    R, rmsdval = align.rotation_matrix(ref0, refu0)
    refCA.translate(-refCA.center_of_mass())
    refCA.rotate(R)
    ref0 = refCA.positions

    CA = u.select_atoms('name CA and %s'%selection)
    noh = u.select_atoms('protein and not type H and %s'%selection)

    rmsds = []
    t = []
    nframe = len(u.trajectory)
    nskip =  int(nframe/nsample)
    
    frame = 0
    for ts in u.trajectory:
        if frame % nskip == 0:
            print(ts.time) 
            traj0 = CA.positions - CA.center_of_mass()
            R, rmsdval = align.rotation_matrix(traj0, ref0)
            noh.translate(-CA.center_of_mass())
            noh.rotate(R)
            noh.write('pdb_noh/%s_65_230/%s_frame%d.pdb'%(mutant,mutant,frame))
        frame += 1
'''
r = []
for mutant in mutants:
    for i in range(nsample+1):
        frame = i*100
        trajwhole = mda.Universe('pdb_noh/%s_65_230/%s_frame%d.pdb'%(mutant,mutant,frame)).select_atoms('protein and not type H')
        comp = ComplexFiltration()
        comp.vertices = trajwhole.positions
        r.append(np.amax(np.linalg.norm(comp.vertices,axis=1)))
rmax = np.amax(r)
print('Rmax = %.3f'%rmax)

for mutant in mutants:
    for i in range(nsample+1):
        frame = i*100
        trajwhole = mda.Universe('pdb_noh/%s_65_230/%s_frame%d.pdb'%(mutant,mutant,frame)).select_atoms('protein and not type H')
        sys.stdout.write('Constructing topology for %s for Frame %d...\r'%(mutant,frame))
        sys.stdout.flush()
        comp = ComplexFiltration()
        comp.vertices = trajwhole.positions
        comp.calc_distance_matrix()
        edges, distances = comp.get_edge_list(radius=4)
        faces = comp.edge_to_face_list(edges=edges)
        comp.vertices /= rmax
        comp.write_off_file(edges=edges,faces=faces,filename='off_noh/%s_65_230_4/%s_4_frame%d.off'%(mutant,mutant,frame))
    

