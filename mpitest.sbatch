#!/bin/tcsh
#
#SBATCH -p normal
#
# S B A T C H -N 1
# S B A T C H --exclusive
#SBATCH -n 16
#SBATCH --ntasks-per-node 16
#SBATCH --output=/home/nsloat/swellix/mpitest/42mer_bundle_%J.o
#SBATCH -t 1200
#SBATCH --job-name=bundle_42mer
#SBATCH --mail-user=nsloat@ou.edu
#
module load mpi/openmpi/1.10.1/intel
#
#date
/usr/bin/time -p mpirun -np 8 /home/nsloat/swellix/swellix-mpi.exe -i /home/nsloat/swellix/42.seq \
								   -M "((x))&CUxAG" -s a -b
#date
