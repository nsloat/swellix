#!/bin/tcsh
#
#SBATCH -p normal
#
# S B A T C H -N 1
#SBATCH --exclusive
#SBATCH -n 1
# S B A T C H --ntasks-per-node 8
#SBATCH --output=/home/nsloat/swellix/perf_exp2_withStats/testNOSTATS.o
#SBATCH -t 0-48:00
#SBATCH --job-name=nostatsComp
#SBATCH --mail-user=nsloat@ou.edu
#
#module load mpi/openmpi/1.10.1/intel
#
#date
/usr/bin/time /home/nsloat/swellix/swellix.exe -i /home/nsloat/swellix/perf_exp2_withStats/data/input.174 -b -l 5
#date
