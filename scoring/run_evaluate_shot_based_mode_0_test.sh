#!/bin/bash

NUMBER_RESULTS=100
GT_FILE=test_dataset_ground_truth_local.txt
RESULTS_SCENES_FILE=../retriever/test/results/SCFV_shots/gaussians_512/ws_mode_0/ws_thresh_8/min_num_words_visited_0/shot_mode_0/fpshot_5/shot_thresh_0.8/out_scene_results.txt

# Compose command line
cmd=$(echo \
    python evaluate_scene_retrieval.py \
    $RESULTS_SCENES_FILE \
    $GT_FILE \
    $NUMBER_RESULTS \
    --force)

# Write and execute command
echo $cmd
$cmd
