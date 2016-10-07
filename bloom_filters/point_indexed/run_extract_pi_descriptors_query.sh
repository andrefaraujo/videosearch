#!/bin/bash -e

QUERY_LIST=test_query_list.txt
ls ../../indexer/test_query/*jpg > $QUERY_LIST

# Parameters for frame-based index
GDINDEX_PATH=../../indexer/global_descriptors/trained_parameters
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

# Compose output index name
out_index=${QUERY_LIST%.txt}.${LD_NAME}_fv_point_idx_k$CENTROIDS

# Command line
cmd=$(echo time \
           ./index_dataset \
           -i ${QUERY_LIST} \
           -o $out_index \
           -r $GDINDEX_PATH \
           -c $CENTROIDS \
           -l $LD_MODE \
           -v $VERBOSE)

# Write and execute command
echo $cmd
$cmd

