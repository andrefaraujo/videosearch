#!/bin/bash -e

QUERY_LIST=test_query_list.txt
ls ../../indexer/test_query/*jpg > $QUERY_LIST

# Parameters for frame-based index
CENTROIDS=512
LD_MODE=1
VERBOSE=1

if [ $LD_MODE -eq 0 ]; then
    LD_NAME=sift
elif [ $LD_MODE -eq 1 ]; then
    LD_NAME=siftgeo
else
    echo "Unrecognized LD_NAME"
    exit
fi

in_name=${QUERY_LIST%.txt}.${LD_NAME}_fv_point_idx_k${CENTROIDS}
out_name=${QUERY_LIST%.txt}.${LD_NAME}_bfv_point_idx_k${CENTROIDS}
    
# Command line
cmd=$(echo time \
           ./binarize_index \
           -i $in_name \
           -o $out_name \
           -v $VERBOSE)

# Write and execute command
echo $cmd
$cmd

