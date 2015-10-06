#!/bin/bash

# Get list of keyframe lists
KEYFRAMES_LISTS=test_keyframe_lists.txt
ls ../test_db/*.txt > $KEYFRAMES_LISTS

# Parameters for index
GDINDEX_PATH=trained_parameters
CENTROIDS=512
LD_MODE=0
VERBOSE=1

# Shot parameters
SHOT_MODE=0
SHOT_KEYF=5
SHOT_THRESH=0.8

# Loop over each video and get shot-based index
# for each of them
for list in `cat $KEYFRAMES_LISTS`; do

    # Get shot file
	shot_res_file=${list%.txt}.shot_t$SHOT_THRESH

    # Compose output index name
    out_index=${list%.txt}.sift_scfv_idx_k${CENTROIDS}_shot_t${SHOT_THRESH}_n${SHOT_KEYF}_m${SHOT_MODE}

    # Command line
    cmd=$(echo time \
        ./index_dataset \
	    -i $list \
        -o $out_index \
        -r $GDINDEX_PATH \
        -c $CENTROIDS \
        -l $LD_MODE \
        -v $VERBOSE \
        -m $SHOT_MODE \
        -k $SHOT_KEYF \
        -s $shot_res_file)

    # Write and execute command
    echo $cmd
    $cmd
done
