#!/bin/bash

NUMBER_RESULTS=100
GT_FILE=test_dataset_ground_truth_local.txt
RESULTS_SCENES_FILE=../retriever/test/results/SCFV_scenes/gaussians_scene_512/ws_mode_0/ws_thresh_scene_6/min_num_words_visited_0/scene_mode_1/fpscene_-1/num_scene_rerank_100/gaussians_shot_512/ws_thresh_shot_7/shot_mode_1/fpshot_-1/shot_thresh_0.8/out_scene_results.txt

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
