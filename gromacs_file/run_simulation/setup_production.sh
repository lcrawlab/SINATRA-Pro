#!/bin/bash
module load mpi
module load cuda/7.5.18
source /users/wtang8/program/plumed-2.6.0/sourceme.sh

gmx="/users/wtang8/program/gromacs/gromacs-2018.8/exec_nompi/bin/gmx"

## set proot to the name of your protein
proot=ubq

this=ptwte
ff=amber99sbdisp
s=1
fileroot=${proot}_${ff}

mdp="mdp/trex_ini.mdp"
gro="data/gro/${fileroot}_${this}_eq.part0001.gro"

## set top to the name of the .top file you have generated
top="${proot}_${ff}_tip4pd_ions.top"

$gmx grompp -v -c $gro -o data/tpr/${proot}_${ff}_${this}_s${s}.tpr -f $mdp -p $top -maxwarn 1

rm \#*
