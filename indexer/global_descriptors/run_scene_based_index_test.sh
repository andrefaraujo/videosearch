#!/bin/bash -e

# Get list of keyframe lists
KEYFRAMES_LISTS=test_keyframe_lists.txt
ls ../test_db/*.txt > $KEYFRAMES_LISTS

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

# Loop over each video and get shot-based index
# for each of them
for list in `cat $KEYFRAMES_LISTS`; do

    # Compose output index name
    out_index=${list%.txt}.${LD_NAME}_scfv_idx_k${CENTROIDS}_scene_n${SCENE_KEYF}_m${SCENE_MODE}

    # Command line
    cmd=$(echo time \
        ./index_dataset \
	    -i $list \
        -o $out_index \
        -r $GDINDEX_PATH \
        -c $CENTROIDS \
        -l $LD_MODE \
        -v $VERBOSE \
        -m $SCENE_MODE \
        -k $SCENE_KEYF \
        --single_shot)

    # Write and execute command
    echo $cmd
    $cmd
done
