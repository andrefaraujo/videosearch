#!/bin/bash

# Get list of keyframe lists
KEYFRAMES_LISTS=test_keyframe_lists.txt
ls ../test_db/*.txt > $KEYFRAMES_LISTS

# Parameters for frame-based index
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

# Loop over each video and get frame-based index
# for each of them
for list in `cat $KEYFRAMES_LISTS`; do

    # Compose output index name
    out_index=${list%.txt}.${LD_NAME}_scfv_idx_k$CENTROIDS

    # Command line
    cmd=$(echo time \
        ./index_dataset \
	    -i $list \
        -o $out_index \
        -r $GDINDEX_PATH \
        -c $CENTROIDS \
        -l $LD_MODE \
        -v $VERBOSE)

    # Write and execute command
    echo $cmd
    $cmd
done
