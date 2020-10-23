#!/bin/bash

proot=ubq
ff=amber99sbdisp
this=ptwte

module load mpi
module load cuda/7.5.18
source /users/wtang8/program/plumed-2.6.0/sourceme.sh

gmx="/users/wtang8/program/gromacs/gromacs-2018.8/exec_nompi/bin/gmx"
end=100000

$gmx convert-tpr -s data/tpr/${proot}_${ff}_${this}_s1.tpr -until $end \
           -o data/tpr/${proot}_${ff}_${this}_s1.tpr

