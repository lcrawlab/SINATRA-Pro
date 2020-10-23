#!/bin/bash

module load mpi
module load cuda
source /users/wtang8/program/plumed-2.6.0/sourceme.sh

gmx="/users/wtang8/program/gromacs/gromacs-2018.8/exec_nompi/bin/gmx"

## set proot to the name of your protein
proot=ubq

this=ptwte
ff=amber99sbdisp
s=1

#for((i=0;i<nrep;i++))
#do

mdp="mdp_eq/trex_ini.mdp"

## Set gro to equal the name of your conformation file ##
gro="data/gro/${proot}_${ff}_tip4pd_npt.gro"

## Set top to the name of the .top file you have generated
top="${proot}_${ff}_tip4pd_ions.top"

$gmx grompp -v -c $gro -o data/tpr/${proot}_${ff}_${this}_eq.tpr -f $mdp -p $top -maxwarn 1

#done
rm \#*
