#!/bin/bash -e

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

# Scene parameters
SCENE_MODE=1
SCENE_KEYF=-1

# Get list of indexes
SCENE_INDEXES=test_scene_based_index_lists.txt
ls ../test_db/*.${LD_NAME}_scfv_idx_k${CENTROIDS}_scene_n${SCENE_KEYF}_m${SCENE_MODE} > $SCENE_INDEXES

# Compose output index name
OUT_INDEX=${SCENE_INDEXES%.txt}.${LD_NAME}_scfv_idx_k${CENTROIDS}_scene_n${SCENE_KEYF}_m${SCENE_MODE}

# Command line
cmd=$(echo time \
    ./join_indexes \
	-i $SCENE_INDEXES \
    -o $OUT_INDEX \
    -r $GDINDEX_PATH \
    -c $CENTROIDS \
    -l $LD_MODE \
    -t $THREADS \
    -v $VERBOSE)

# Write and execute command
echo $cmd
$cmd
