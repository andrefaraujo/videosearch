#!/bin/bash

# Get list of queries
IMAGE_LIST=test_query_list.txt
ls ../test_query/*.jpg > $IMAGE_LIST

# Parameters for query index
GDINDEX_PATH=trained_parameters
CENTROIDS=512
LD_MODE=0
VERBOSE=1

# Compose output index name
out_index=${IMAGE_LIST%.txt}.sift_scfv_idx_k$CENTROIDS

# Command line
cmd=$(echo time \
    ./index_dataset \
	-i $IMAGE_LIST \
    -o $out_index \
    -r $GDINDEX_PATH \
    -c $CENTROIDS \
    -l $LD_MODE \
    -v $VERBOSE)

# Write and execute command
echo $cmd
$cmd
