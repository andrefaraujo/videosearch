#!/bin/bash

# Number of threads to use
THREADS=10

# Parameters for index
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

# Shot parameters
SHOT_MODE=0
SHOT_KEYF=5
SHOT_THRESH=0.8

# Get list of indexes
SHOT_INDEXES=test_shot_based_index_lists.txt
ls ../test_db/*.${LD_NAME}_scfv_idx_k${CENTROIDS}_shot_t${SHOT_THRESH}_n${SHOT_KEYF}_m${SHOT_MODE} > $SHOT_INDEXES

# Compose output index name
OUT_INDEX=${SHOT_INDEXES%.txt}.${LD_NAME}_scfv_idx_k${CENTROIDS}_shot_t${SHOT_THRESH}_n${SHOT_KEYF}_m${SHOT_MODE}

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
