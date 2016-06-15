#!/bin/sh

[ -z "${1}" ] && {
  echo "Example usage: ${0} stmv.conf"
  exit
}

#rm -rf $3/structs/$2
#mkdir $3/structs/$2
rm -rf $3/labeled/$2
mkdir $3/labeled/$2
python generate_structures.py $1 $2 $3

# The labeling process used to be run here and with a script in Python.
# It was moved to the label executable and is now run each time a .struct file is produced.
# This allows the overall disk footprint at any given time to be decreased by deleting the .struct
# file as soon as it's processed. In addition, the C implementation of this labeling program is much 
# faster.
#python label.py $1 $2
