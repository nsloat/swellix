#!/bin/tcsh
#
#SBATCH -p normal
#
# S B A T C H -N 1
#SBATCH --exclusive
#SBATCH -n 112
#SBATCH --ntasks-per-node 16
#SBATCH --output=/scratch/nsloat/mpitest/par_HERV244_l3_%J.o
#SBATCH -t 2-00:00
#SBATCH --job-name=HERV244_par
#SBATCH --mail-user=nsloat@ou.edu
#
module load mpi/openmpi/1.10.1/intel
#
date
/usr/bin/time -p mpirun -np 112 /home/nsloat/swellix/dev-swellix-mpi.exe \
	-i /home/nsloat/seq/HERV_244.seq -b -l 3 -d 1
date
