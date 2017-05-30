#!/bin/bash

source ~/.bashrc

EXPDIR="/u/eot/sloat/swellix/perf_exp2_mpi_noStats"

S="/u/eot/sloat/swellix/dev-swellix-mpi.exe"
S_SBATCH=$EXPDIR"/mpiexp.pbs"

#s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;
	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.14 -b -d 1 -l "$i
    NODES=1
    NUMPROC=7
    PPN=7
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/14/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
    qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.28 -b -l "$i
    NUMPROC=28
    NODES=1
    PPN=28
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/28/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.42 -b -l "$i
    NUMPROC=42
    NODES=2
    PPN=21
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/42/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.50 -b -l "$i
    NUMPROC=50
    NODES=2
    PPN=25
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/50/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.71 -b -l "$i
    NUMPROC=71
    NODES=3
    PPN=25
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/71/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.74 -b -l "$i
    NUMPROC=74
    NODES=3
    PPN=25
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/74/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.76 -b -l "$i
    NUMPROC=76
    NODES=3
    PPN=26
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/76/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.161 -b -l "$i
    NUMPROC=161
    NODES=6
    PPN=27
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/161/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.167 -b -l "$i
    NUMPROC=167
    NODES=6
    PPN=28
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/167/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.171 -b -l "$i
    NUMPROC=171
    NODES=6
    PPN=29
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/171/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.174 -b -l "$i
    NUMPROC=174
    NODES=6
    PPN=29
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/174/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.197 -b -l "$i
    NUMPROC=197
    NODES=7
    PPN=29
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/197/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done

	for i in $(seq 2 10)
	do

		#Configure and run Swellix for sequence i
		ARGS="-i "$EXPDIR"/data/input.247 -b -l "$i
    NUMPROC=247
    NODES=8
    PPN=31
		sed 's:\*NODES\*:'$NODES':g;s:\*NP\*:'$NUMPROC':g;s:\*NPERNODE\*:'$PPN':g;s:\*OUTPATH\*:'$EXPDIR"/data/247/length"$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
		
		qsub swellix.sbatch

	done
rm swellix.sbatch
