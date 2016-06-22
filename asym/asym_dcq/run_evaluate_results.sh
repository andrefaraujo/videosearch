#!/bin/bash

NUMBER_RESULTS=100
GT_FILE=ground_truth_public.txt
RESULTS_FILE=results/C_0005/SCFV_images/gaussians_512/ws_mode_1/ws_thresh_0.0004/min_num_words_visited_0/out_results.txt

# Compose command line
cmd=$(echo \
    python ../../scoring/evaluate_scene_retrieval.py \
    $RESULTS_FILE \
    $GT_FILE \
    $NUMBER_RESULTS \
    --force)

# Write and execute command
echo $cmd
$cmd
