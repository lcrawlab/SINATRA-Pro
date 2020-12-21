#!/bin/python3
import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import RMSF

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib as mpl

mpl.rcParams['font.sans-serif'] = 'Arial'
mpl.rcParams['font.family'] = 'sans-serif'

directory='rmsf/'
mutant = "R164G"

u = mda.Universe('%s/md_0_1.gro'%mutant,'%s/md_0_1_noPBC.xtc'%mutant)
protein = u.select_atoms('protein and resid 65:213')

# 1) need a step to center and make whole: this trajectory
#    contains the protein being split across periodic boundaries
#
# TODO

# 2) fit to the initial frame to get a better average structure
#    (the trajectory is changed in memory)
prealigner = align.AlignTraj(u, u, select="protein and name CA and resid 65:213", in_memory=True).run()

# 3) reference = average structure
reference_coordinates = u.trajectory.timeseries(asel=protein).mean(axis=1)
# make a reference structure (need to reshape into a 1-frame "trajectory")
reference = mda.Merge(protein).load_new(reference_coordinates[:, None, :], order="afc")

aligner = align.AlignTraj(u, reference, select="protein and name CA and resid 65:213", in_memory=True).run()

calphas = protein.select_atoms("name CA and resid 65:213")
rmsfer = RMSF(calphas, verbose=True).run()

np.savetxt(directory+'%s_resnums.dat'%mutant,calphas.resnums,fmt='%d')
np.savetxt(directory+'%s_rmsf.dat'%mutant,rmsfer.rmsf,fmt='%.6f')


