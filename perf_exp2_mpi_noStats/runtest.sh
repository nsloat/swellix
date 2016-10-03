#!/bin/bash

source ~/.bashrc

EXPDIR="/home/nsloat/swellix/perf_exp2_mpi_noStats"

S="/home/nsloat/swellix/dev-swellix-mpi.exe"
S_SBATCH=$EXPDIR"/mpiexp_template.sbatch"

#s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;
	for i in $(seq 2 2)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.14 -b -d 1 -l "$i
    NUMPROC=7
    PPN=7
		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/14/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		sbatch swellix.sbatch

	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.28 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/28/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.42 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/42/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.50 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/50/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.71 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/71/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.74 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/74/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.76 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/76/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.161 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/161/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.167 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/167/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.171 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/171/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.174 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/174/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.197 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/197/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
#
#	for i in $(seq 2 10)
#	do
#
#		#Configure and run Swellix for sequence i
#		ARGS="-i "$EXPDIR"/data/input.247 -b -l "$i
#    NUMPROC=16
#    PPN=16
#		sed 's:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/247/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#		
#		sbatch swellix.sbatch
#
#	done
rm swellix.sbatch
