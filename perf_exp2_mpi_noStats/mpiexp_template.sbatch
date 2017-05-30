#!/bin/tcsh
#
#SBATCH -p normal
#
# S B A T C H -N 1
#SBATCH --exclusive
#SBATCH -n *NP*
#SBATCH --ntasks-per-node *NPERNODE*
#SBATCH --output=*OUTPATH*
#SBATCH -t 0-48:00
#SBATCH --job-name=pe2_mpi
#SBATCH --mail-user=nsloat@ou.edu
#
module load mpi/openmpi/1.10.1
#
#date
/usr/bin/time -p mpirun -np *NP* *EXEC* *ARGS*
#date
