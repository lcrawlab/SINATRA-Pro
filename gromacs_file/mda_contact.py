#!/bin/python
import sys
import numpy as np
import MDAnalysis as mda
from MDAnalysis.lib.nsgrid import FastNS, NSResults
from MDAnalysis.analysis import contacts, distances

#trajfile = 'prot.xtc'
topfile = sys.argv[1] #'prot.gro'
trajfile = sys.argv[2] #'data/xtc_prot/FUS_120-163_FUS_454-501_amber99sbwsSTQ_ptwte_s1_nd0_eq.xtc'

u = mda.Universe(topfile,trajfile)
proteinNOH = u.select_atoms('protein and not type H and resid 65:213')
atomind = proteinNOH.atoms.ix
resid = np.zeros(len(proteinNOH),dtype=int)
residues = proteinNOH.residues
resnum = 0
for a in residues:
    for ind in a.atoms.ix:
        arr = np.where(atomind == ind)
        if len(arr) == 1:
            resid[arr[0]] = resnum
            print(resnum,a.ix)
    resnum += 1
print(residues)
exit()
frame = 0
nframe = len(u.trajectory)
matrix = np.zeros((resnum,resnum,nframe),dtype=bool)

for ts in u.trajectory:
    sys.stdout.write('\rFrame %d'%frame)
    sys.stdout.flush()
    nsr = FastNS(cutoff=6.0,coords=proteinNOH.positions,box=u.dimensions,pbc=True,max_gridsize=20000)
    result = nsr.self_search()
    nb = result.get_pairs()[::2]
    for pair in nb:
        matrix[resid[pair[0]],resid[pair[1]],frame] = True
    frame += 1

np.save(sys.argv[3],matrix)

