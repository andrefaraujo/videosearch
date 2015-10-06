#!/bin/bash

NUMBER_RESULTS=100
RESULTS_FILE=../retriever/test/results/SCFV_shots/gaussians_512/ws_mode_0/ws_thresh_7/min_num_words_visited_0/shot_mode_1/fpshot_-1/shot_thresh_0.8/out_results.txt
RESULTS_SCENES_FILE=../retriever/test/results/SCFV_shots/gaussians_512/ws_mode_0/ws_thresh_7/min_num_words_visited_0/shot_mode_1/fpshot_-1/shot_thresh_0.8/out_scene_results.txt

# Compose command line
cmd=$(echo \
    python convertKeyframeResultsToSceneRetrievalResults.py \
    $RESULTS_FILE \
    $RESULTS_SCENES_FILE \
    $NUMBER_RESULTS)

# Write and execute command
echo $cmd
$cmd
