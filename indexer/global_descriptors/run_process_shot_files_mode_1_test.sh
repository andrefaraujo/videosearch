#!/bin/bash -e

# Get list of keyframe lists
KEYFRAMES_LISTS=test_keyframe_lists.txt
ls ../test_db/*.txt > $KEYFRAMES_LISTS

# Parameters 
CENTROIDS=512
SHOT_MODE=1
FPSHOT=-1
SHOT_THRESH=0.8

OUT_NAME=test_shot_mode_1_first_frames.txt

# Command line
cmd=$(echo python \
    ./generate_shot_based_aux_files.py \
    $KEYFRAMES_LISTS \
    $CENTROIDS \
    $SHOT_THRESH \
    $FPSHOT \
    $SHOT_MODE \
    dummy \
    $OUT_NAME)

# Write and execute command
echo $cmd
$cmd
