#!/bin/bash

# Number of threads to use
THREADS=10

# Parameters for index
GDINDEX_PATH=trained_parameters
CENTROIDS=512
LD_MODE=0
VERBOSE=1

# Shot parameters
SHOT_MODE=1
SHOT_KEYF=-1
SHOT_THRESH=0.8

# Get list of indexes
SHOT_INDEXES=test_shot_based_index_lists.txt
ls ../test_db/*.sift_scfv_idx_k${CENTROIDS}_shot_t${SHOT_THRESH}_n${SHOT_KEYF}_m${SHOT_MODE} > $SHOT_INDEXES

# Compose output index name
OUT_INDEX=${SHOT_INDEXES%.txt}.sift_scfv_idx_k${CENTROIDS}_shot_t${SHOT_THRESH}_n${SHOT_KEYF}_m${SHOT_MODE}

# Command line
cmd=$(echo time \
    ./join_indexes \
	-i $SHOT_INDEXES \
    -o $OUT_INDEX \
    -r $GDINDEX_PATH \
    -c $CENTROIDS \
    -l $LD_MODE \
    -t $THREADS \
    -v $VERBOSE)

# Write and execute command
echo $cmd
$cmd
