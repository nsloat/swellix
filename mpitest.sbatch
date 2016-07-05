#!/bin/tcsh
#
#SBATCH -p normal
#
# S B A T C H -N 1
# S B A T C H --exclusive
#SBATCH -n 8
#SBATCH --ntasks-per-node=8
#SBATCH --output=/home/nsloat/swellix/swellix/mpitest/test_%J.o
#SBATCH -t 1200
#SBATCH --job-name=slurmMPItest
#SBATCH --mail-user=nsloat@ou.edu
#
module load mpi/openmpi/1.10.1/intel
#
#date
/usr/bin/time -p mpirun -np 8 /home/nsloat/swellix/swellix/swellix-mpi.exe -i /home/nsloat/swellix/swellix/14.seq -b 
#date
