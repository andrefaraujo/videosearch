#!/bin/bash -e

# Get list of videos
VIDEO_LIST=test_video_list.txt
ls ../test_db/test_video_*.mp4 > $VIDEO_LIST

# Parameters 
SHOT_THRESH=0.8

OUT_NAME=test_group_lists_rerank.txt

# Command line
cmd=$(echo python \
    ./generate_scene_based_rerank_aux_files.py \
    $VIDEO_LIST \
    $SHOT_THRESH \
    $OUT_NAME)

# Write and execute command
echo $cmd
$cmd
