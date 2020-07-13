#!/bin/bash
#############
# Arguments #
#############

cmp=$1
it=$2
job=$3
wall_time=$4

#########
# Start #
#########

#echo $cmp
file=sub_${it}/${cmp}_${job}_sub.sh

bsub -W ${wall_time} -n 1 -J ${cmp}_${it}_${job} -e out_${it}/${cmp}_${job}.err -o  out_${it}/${cmp}_${job}.out < $file

ret=$?
if [ $ret -ne 0 ]; then
    exit 2
fi

