#!/bin/bash

source ~/.bashrc

EXPDIR=$(pwd)

S="/u/eot/sloat/swellix/swellix.exe"
S_PBS="/u/eot/sloat/swellix/template.pbs"

for DIR in $(ls -d ./*/)
  do
    INFILE=$(ls $DIR | grep ".seq")
    WALLTIME="48:00:00"
    OUTPATH=$EXPDIR"/"$DIR

#    NUMPROC=$(($(wc -m $EXPDIR"/"$DIR"/"$INFILE | awk '{print $1}') - 1))
    NUMPROC=1
    PPN=1
    NODES=1
    for MHL in $(seq 2 6)
      do
      SWELLARGS="-b -l "$MHL" -i "$EXPDIR"/"$DIR$INFILE
      JOBNAME="px4_mhl"$MHL"_"$(echo $INFILE | awk -F "." '{ print $1 }')
#    CMS=$(echo $SEQDIR | awk '{print gensub(/^([0-9]*)_CMs.*$/, "\\1", "g")}')
#               echo $SEQDIR, $CMS
#    i=$(ls -l "./data/"$SEQDIR | grep .swlx | wc -l)
#               i=$(ls -l "./data/5_CMs" | grep .swlx | wc -l)
#               CMS=5

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
      qsub swellix.pbs

    done
  done

rm swellix.pbs

