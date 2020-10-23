#!/bin/bash -l
#----------------------------------------------------
# Example SLURM job script to run MPI applications on 
# OSCAR system.
#----------------------------------------------------

# Request an hour of runtime:
#SBATCH --time=24:00:00

# Use 2 nodes with 8 tasks each, for 16 MPI tasks:
#SBATCH --nodes=2
#SBATCH --tasks-per-node=16
#SBATCH --mem=32G
#SBATCH --begin=now+4hour

# Specify a job name:
#SBATCH -J ubq_WT_run

# Specify an output file
#SBATCH -o MyMPIJob-%j.out
#SBATCH -e MyMPIJob-%j.out

module load mpi
module load cuda
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

options="-maxh 24 -npme 0 -ntomp 1 \
-v -s data/tpr/${fileroot}_${this}_s1.tpr \
-x data/xtc/${fileroot}_${this}_s${s}.xtc \
-o data/trr/${fileroot}_${this}_s${s}.trr \
-c data/gro/${fileroot}_${this}_s${s}.gro \
-e data/edr/${fileroot}_${this}_s${s}.edr \
-g data/log/${fileroot}_${this}_s${s}.log \
-cpo data/cpt/${fileroot}_${this}_s${s}.cpt \
-cpi data/cpt/${fileroot}_${this}_s${s}.cpt -noappend"

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

srun --mpi=pmi2 $application $options > outfile_${proot} 2>&1

echo Time is `date`
