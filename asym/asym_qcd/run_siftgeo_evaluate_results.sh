#!/bin/bash -e

NUMBER_RESULTS=100
GT_FILE=ground_truth_public.txt
RESULTS_FILE=results/C_0005/FV_shots/gaussians_512/ws_mode_1/ws_thresh_0.0008/min_num_words_visited_0/asym_scoring_mode_1/out_results.txt

# Compose command line
cmd=$(echo \
    python ../../scoring/evaluate_scene_retrieval.py \
    $RESULTS_FILE \
    $GT_FILE \
    $NUMBER_RESULTS)

# Write and execute command
echo $cmd
$cmd
