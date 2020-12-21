#!/bin/python3

import numpy as np
import MDAnalysis as mda

mutant = "pD58"

u = mda.Universe('%s/md_0_1.gro'%mutant,'%s/md_0_1_noPBC.xtc'%mutant)
u.trajectory[-1]
protein = u.select_atoms('protein')
protein.write("%s_last.gro"%mutant)

