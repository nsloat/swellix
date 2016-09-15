#!/bin/tcsh
#
#SBATCH -p normal
#
# S B A T C H -N 1
#SBATCH --exclusive
#SBATCH -n 2
#SBATCH --ntasks-per-node 2
#SBATCH --output=/home/nsloat/swellix/mpitest/bundMPI_%J.o
#SBATCH -t 1200
#SBATCH --job-name=bundMPItest
#SBATCH --mail-user=nsloat@ou.edu
#
module load mpi/openmpi/1.10.1/intel
#
#date
/usr/bin/time -p mpirun -np 2 /home/nsloat/swellix/dev-swellix-mpi.exe -i /home/nsloat/swellix/28.seq -b
#date
