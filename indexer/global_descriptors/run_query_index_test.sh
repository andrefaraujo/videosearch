#!/bin/bash

# Get list of queries
IMAGE_LIST=test_query_list.txt
ls ../test_query/*.jpg > $IMAGE_LIST

# Parameters for query index
GDINDEX_PATH=trained_parameters
CENTROIDS=512
LD_MODE=0
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
out_index=${IMAGE_LIST%.txt}.${LD_NAME}_scfv_idx_k$CENTROIDS

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
