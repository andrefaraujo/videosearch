#!/bin/bash

# Number of threads to use
THREADS=10

# Parameters for frame-based index
GDINDEX_PATH=trained_parameters
CENTROIDS=512
LD_MODE=0
VERBOSE=1

# Get list of indexes
KEYFRAMES_INDEXES=test_frame_based_index_lists.txt
ls ../test_db/*.sift_scfv_idx_k$CENTROIDS > $KEYFRAMES_INDEXES

# Compose output index name
OUT_INDEX=${KEYFRAMES_INDEXES%.txt}.sift_scfv_idx_k$CENTROIDS

# Command line
cmd=$(echo time \
    ./join_indexes \
	-i $KEYFRAMES_INDEXES \
    -o $OUT_INDEX \
    -r $GDINDEX_PATH \
    -c $CENTROIDS \
    -l $LD_MODE \
    -t $THREADS \
    -v $VERBOSE)

# Write and execute command
echo $cmd
$cmd
