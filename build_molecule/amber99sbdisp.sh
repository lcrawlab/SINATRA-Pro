#!/bin/bash -l
#----------------------------------------------------
# Example SLURM job script to run MPI applications on 
# OSCAR system.
#----------------------------------------------------

# Request an hour of runtime:
#SBATCH --time=00:20:00

# Use 2 nodes with 8 tasks each, for 16 MPI tasks:
#SBATCH --nodes=2
#SBATCH --tasks-per-node=6
#SBATCH --mem=8G

# Specify a job name:
#SBATCH -J 1UBQ_WT_build

# Specify an output file
#SBATCH -o MyMPIJob-%j.out
#SBATCH -e MyMPIJob-%j.out

module load mpi
module load cuda
source /users/wtang8/program/plumed-2.6.0/sourceme.sh
gromacs="/users/wtang8/program/gromacs/gromacs-2018.8/exec_plumed/bin/gmx_mpi"
gmxnompi="/users/wtang8/program/gromacs/gromacs-2018.8/exec_nompi/bin/gmx"

proot=ubq

ff="amber99sbdisp"
water="tip4pd"
waterG="tip4p"
fileroot="${proot}_${ff}_${water}"
#pdbfile="data/pdb/${proot}.gro"
pdbfile="1ubq_prot.pdb"
boxlen=6 # nm

PDB2GMX="$gmxnompi pdb2gmx"
EDITCONF="$gmxnompi editconf"
GROMPP="$gmxnompi grompp"
MDRUN="$gromacs mdrun"
GENBOX="$gmxnompi genbox"
GENION="$gmxnompi genion"
SOLVATE="$gmxnompi solvate"

printf '1\n1\n' | $PDB2GMX -f ${pdbfile} \
	-p ${fileroot}.top \
	-o data/gro/${fileroot}_vac.gro \
	-ignh -water select \
	-i ${fileroot}_posre.itp -renum -ss -ter

## max protein principal axis ~26.0 A
$EDITCONF -bt octahedron -f data/gro/${fileroot}_vac.gro \
	-o data/gro/${fileroot}_ed.gro \
	-box $boxlen

$GROMPP -v -f mdp/mini.mdp -c data/gro/${fileroot}_ed.gro \
        -o data/tpr/${fileroot}_vac_min.tpr -p ${fileroot}.top

srun --mpi=pmi2 $MDRUN -v -npme 0 \
        -s data/tpr/${fileroot}_vac_min.tpr \
        -o data/trr/${fileroot}_vac_min.trr \
        -c data/gro/${fileroot}_vac_min.gro \
        -e data/edr/${fileroot}_vac_min.edr \
        -g data/output/${fileroot}_vac_min.log

$SOLVATE -cp data/gro/${fileroot}_vac_min.gro \
	-cs ${waterG}.gro -o data/gro/${fileroot}_wat.gro \
	-p ${fileroot}.top
cp ${fileroot}.top backup.top

sed -i 's/ MW/MW4/g' data/gro/${fileroot}_wat.gro
sed -i 's/HW2/HW3/g' data/gro/${fileroot}_wat.gro
sed -i 's/HW1/HW2/g' data/gro/${fileroot}_wat.gro

$GROMPP -v -f mdp/mini.mdp -c data/gro/${fileroot}_wat.gro \
	-o data/tpr/${fileroot}_mini.tpr -p ${fileroot}.top

srun --mpi=pmi2 $MDRUN -v -ntomp 1 -npme 0 \
        -s data/tpr/${fileroot}_mini.tpr \
	-o data/trr/${fileroot}_mini.trr \
	-c data/gro/${fileroot}_wat_min.gro \
	-e data/edr/${fileroot}_mini.edr \
	-g data/output/${fileroot}_mini.log

$GROMPP -v -f mdp/pme.mdp -c data/gro/${fileroot}_wat_min.gro \
         -o pme.tpr -p ${fileroot}.top -maxwarn 2
cp ${fileroot}.top ${fileroot}_ions.top

## Add salt
echo 'SOL' | $GENION -s pme.tpr -p ${fileroot}_ions.top -o data/gro/${fileroot}_ions.gro -neutral -conc 0.1 
nna=`grep 'NA' ${fileroot}_ions.top | tail -n 1 | grep -o '[0-9]*'`
ncl=`grep 'CL' ${fileroot}_ions.top | tail -n 1 | grep -o '[0-9]*'`
echo $nna $ncl
excessions=$[$[$nna-$ncl]/2]
excessions=${excessions#-}
nna=$[$nna-$excessions]
ncl=$[$ncl-$excessions]
echo $nna $ncl
cp ${fileroot}.top ${fileroot}_ions.top
echo 'SOL' | $GENION -s pme.tpr -p ${fileroot}_ions.top -o data/gro/${fileroot}_ions.gro -np $nna -nn $ncl 

$GROMPP -v -f mdp/nvt.mdp -c data/gro/${fileroot}_ions.gro \
        -o data/tpr/${fileroot}_nvt.tpr -p ${fileroot}_ions.top

srun --mpi=pmi2 $MDRUN -v -ntomp 1 -npme 0 -ntomp_pme 1\
        -s data/tpr/${fileroot}_nvt.tpr \
	-o data/trr/${fileroot}_nvt.trr \
	-c data/gro/${fileroot}_nvt.gro \
	-e data/edr/${fileroot}_nvt.edr \
	-g data/output/${fileroot}_nvt.log

$GROMPP -v -maxwarn 1 -f mdp/npt.mdp -c data/gro/${fileroot}_nvt.gro \
	-o data/tpr/${fileroot}_npt.tpr -p ${fileroot}_ions.top

srun --mpi=pmi2 $MDRUN -v -ntomp 1 -npme 0 -ntomp_pme 1\
        -s data/tpr/${fileroot}_npt.tpr \
        -o data/trr/${fileroot}_npt.trr \
        -c data/gro/${fileroot}_npt.gro \
        -e data/edr/${fileroot}_npt.edr \
 	-x data/xtc/${fileroot}_npt.xtc \
	-cpo data/cpt/${fileroot}_npt.cpt \
        -g data/output/${fileroot}_npt.log

rm data/*/\#*
rm \#*


