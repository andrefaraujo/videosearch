#!/bin/bash -e

# Feature mode: 0=SIFT; 1=SIFTGeo
FEAT_MODE=0
if [ $FEAT_MODE -eq 0 ]; then
    FEAT_NAME=sift
elif [ $FEAT_MODE -eq 1 ]; then
    FEAT_NAME=siftgeo
else
    echo "Unrecognized LD_NAME"
    exit
fi

NUMBER_RESULTS=100
GT_FILE=test_dataset_ground_truth_local.txt
RESULTS_SCENES_FILE=../retriever/test_query_index/results/SCFV_frames/feat_name_${FEAT_NAME}/gaussians_512/ws_mode_0/ws_thresh_10/min_num_words_visited_0/out_scene_results.txt

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
