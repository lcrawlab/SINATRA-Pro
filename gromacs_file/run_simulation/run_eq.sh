#!/bin/bash -l
#----------------------------------------------------
# Example SLURM job script to run MPI applications on 
# OSCAR system.
#----------------------------------------------------

# Request an hour of runtime:
#SBATCH --time=0:30:00

# Use 2 nodes with 8 tasks each, for 16 MPI tasks:
#SBATCH --nodes=2
#SBATCH --tasks-per-node=16

# Specify a job name:
#SBATCH -J ubq_WT_eq

# Specify an output file
#SBATCH -o MyMPIJob-%j.out
#SBATCH -e MyMPIJob-%j.out

module load mpi
module load cuda/7.5.18
source /users/wtang8/program/plumed-2.6.0/sourceme.sh

# Which set?
s=1
# Full path to application + application name
application="/users/wtang8/program/gromacs/gromacs-2018.8/exec/bin/gmx_mpi mdrun"
# Define variables related to protein and ff
proot="ubq"
ff="amber99sbdisp"
fileroot="${proot}_${ff}"
this="ptwte"

options="-maxh 0.5 -npme 0 -ntomp 1 \
-v -s data/tpr/${fileroot}_${this}_eq.tpr \
-x data/xtc/${fileroot}_${this}_eq.xtc \
-o data/trr/${fileroot}_${this}_eq.trr \
-c data/gro/${fileroot}_${this}_eq.gro \
-e data/edr/${fileroot}_${this}_eq.ene \
-g data/log/${fileroot}_${this}_eq.log \
-cpo data/cpt/${fileroot}_${this}_eq.cpt \
-cpi data/cpt/${fileroot}_${this}_eq.cpt -noappend"

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

# Launch the MPI executable named "a.out"

srun --mpi=pmi2 $application $options > outfile_${proot} 2>&1

echo Time is `date`
