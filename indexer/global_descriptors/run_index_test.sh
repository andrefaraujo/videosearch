#!/bin/bash

IMAGE_LIST=../test_all_db_keyframes.txt
OUT_INDEX=../test_all_db_keyframes.sift_scfv_idx_k512
GDINDEX_PATH=trained_parameters
CENTROIDS=512
LD_MODE=0
VERBOSE=1

cmd=$(echo time\
	./index_dataset \
	-i $IMAGE_LIST \
    -o $OUT_INDEX \
    -r $GDINDEX_PATH \
    -c $CENTROIDS \
    -l $LD_MODE \
    -v $VERBOSE)

# Write and execute command
echo $cmd
$cmd
