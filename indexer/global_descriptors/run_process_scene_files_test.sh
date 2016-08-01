#!/bin/bash -e

# Get list of keyframe lists
KEYFRAMES_LISTS=test_keyframe_lists.txt
ls ../test_db/*.txt > $KEYFRAMES_LISTS

# Parameters 
CENTROIDS=512
SCENE_MODE=1
FPSHOT=-1

OUT_NAME=test_scene_mode_1_first_frames.txt

# Command line
cmd=$(echo python \
    ./generate_scene_based_aux_files.py \
    $KEYFRAMES_LISTS \
    $CENTROIDS \
    $FPSHOT \
    $SCENE_MODE \
    dummy \
    $OUT_NAME)

# Write and execute command
echo $cmd
$cmd
