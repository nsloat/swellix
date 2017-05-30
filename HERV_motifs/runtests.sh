#!/bin/bash

source ~/.bashrc

EXPDIR=$(pwd)

S="/u/eot/sloat/swellix/dev-swellix-mpi.exe"
S_PBS="/u/eot/sloat/swellix/template.pbs"

for DIR in $(find . -type d \! -regex '.*HERVS.*' -a \! -regex '\.')
  do
    for MFILE in $(ls $EXPDIR"/"$DIR | grep ".motif")
      do
        MOTIF=$(cat $EXPDIR"/"$DIR"/"$MFILE | awk '{print gensub("&", "\\\\&", "g")}')
        for sequence in $(ls "HERVS" | grep ".seq")
          do
            INFILE=$sequence
            WALLTIME="48:00:00"
            OUTPATH=$EXPDIR"/"$DIR

#    NUMPROC=$(($(wc -m $EXPDIR"/"$DIR"/"$INFILE | awk '{print $1}') - 1))
            MHL=3
            SWELLARGS="-d 1 -b -l "$MHL" -i "$EXPDIR"/HERVS/"$INFILE" -motif "$MOTIF
            JOBNAME=$(echo $MFILE | awk -F "." '{print $1}')"_"$(echo $INFILE | awk -F "." '{ print $1 }')"_dyn"
            SEQLEN=$(echo $JOBNAME | awk -F "_" '{print $3}')
            NUMPROC=$SEQLEN
            PPN=16
            NODES=$((($NUMPROC + 16 - 1)/16))
echo $JOBNAME
echo $NUMPROC
echo $NODES

      #Configure and run Swellix for sequence i
#      ARGS="-k -d 1 -b -l 2 -s a"
            sed 's!\*NODES\*!'"$NODES"'!g;
                s!\*WALLTIME\*!'"$WALLTIME"'!g;
                s!\*JOBNAME\*!'"$JOBNAME"'!g;
                s!\*OUTPATH\*!'"$OUTPATH"'!g;
                s!\*NUMPROC\*!'"$NUMPROC"'!g;
                s!\*PPN\*!'"$PPN"'!g;
                s!\*EXEC\*!'"$S"'!g;
                s!\*SWELLARGS\*!'"$SWELLARGS"'!g' $S_PBS > swellix.pbs
#echo $SWELLARGS
           qsub swellix.pbs
        done  
    done
done
  
rm swellix.pbs

