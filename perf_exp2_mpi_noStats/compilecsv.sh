#!/bin/bash


if [ -f results.csv ]
	then rm results.csv
fi


for j in $(ls -d data/*/ | sed 's:data/\([0-9]\+\)/$:\1:g' | sort -n)
do

	echo "Sequence: "$(cat "data/input."$j) >> results.csv
	echo "MinHelixLen,Bundling(on/off),SequenceLength,# of Bundles,# of Components,# of Structures,Runtime" >> results.csv
	
#	for len in $(seq 2 10)
#	do
		for i in $(ls "data/$j/" | grep ".o" | sort -n)
		do 
			cat "data/$j/$i" > temp.txt
			RUNTIME=$(grep real temp.txt | sed 's/^[^0-9]*\([0-9]\+\.[0-9]\+\).*$/\1/g')
			SEQLEN=$j
			HELIX=$(echo $i | sed 's/^length\([0-9]\+\).*$/\1/g')
			BUNDLING=$(if [ $(grep \ -b\  "data/"$j"/"$i | wc -l) -gt 0 ]; then echo "ON"; else echo "OFF"; fi)
			BUNDLES=$(grep Bundles: temp.txt | sed 's/^[^0-9]\+\([0-9]\+\).*$/\1/g')
			COMPONENTS=$(grep Components: temp.txt | sed 's/^[^0-9]\+\([0-9]\+\).*$/\1/g')
			STRUCTURES=$(grep Bundled\ structures: temp.txt | sed 's/^[^0-9]\+\([0-9]\+\).*$/\1/g')
			echo $HELIX,$BUNDLING,$SEQLEN,$BUNDLES,$COMPONENTS,$STRUCTURES,$RUNTIME >> results.csv
		done
	#done

	echo "" >> results.csv
done

rm temp.txt
