#!/bin/bash

source ~/.bashrc

EXPDIR=$(pwd)

S="/u/eot/sloat/swellix.exe"
S_PBS="/u/eot/sloat/swellix/template.pbs"

        for SEQDIR in $(ls ./data/ | grep _CMs)
        do
                CMS=$(echo $SEQDIR | awk '{print gensub(/^([0-9]*)_CMs.*$/, "\\1", "g")}')
#               echo $SEQDIR, $CMS
                i=$(ls -l "./data/"$SEQDIR | grep .swlx | wc -l)
#               i=$(ls -l "./data/5_CMs" | grep .swlx | wc -l)
#               CMS=5
#  SEQDIR="1_CMs_3"
#  i=2
#  CMS=1
                        #Configure and run Swellix for sequence i
                        ARGS="-k -d 1 -b -l 2 -s a"
                        sed 's:\*SEQDIR\*:'$SEQDIR':g;s:\*MAXMOD\*:'$i':g;s:\*EXECWORKINGDIR\*:/home/nsloat/    swellix/:g;s:\*EXEC\*:'"$S"':g;s:\*ARGS\*:'"$ARGS"':g' $S_SBATCH > swellix.sbatch
#cat swellix.sbatch
                        sbatch swellix.sbatch

        done

rm swellix.sbatch
~

