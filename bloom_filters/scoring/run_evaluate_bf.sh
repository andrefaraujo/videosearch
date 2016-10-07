#!/bin/bash -e

# Number of Gaussians to use in point-indexed
GAUSSIANS=512
# TF-IDF power normalization
ALPHA=0.5
# Feature mode: 0=SIFT; 1=SIFTGeo
FEAT_MODE=1
if [ $FEAT_MODE -eq 0 ]; then
    FEAT_NAME=sift
elif [ $FEAT_MODE -eq 1 ]; then
    FEAT_NAME=siftgeo
else
    echo "Unrecognized LD_NAME"
    exit
fi

NUMBER_RESULTS=100
GT_FILE=../../scoring/test_dataset_ground_truth_local.txt

PATH_PREFIX=test_bf_results
METHOD_PATH_RESULTS=${METHOD_PATH}"/feat_name_"$FEAT_NAME
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/gaussians_"$GAUSSIANS
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/alpha_"$ALPHA
RESULTS_SCENES_FILE=../retriever/${PATH_PREFIX}/results/${METHOD_PATH_RESULTS}/out_results.txt
RESULTS_SCENES_FILE_MODIF=${RESULTS_SCENES_FILE}.modif
sed 's,../../,../,g' $RESULTS_SCENES_FILE > $RESULTS_SCENES_FILE_MODIF

# Compose command line
cmd=$(echo \
          python ../../scoring/evaluate_scene_retrieval.py \
          $RESULTS_SCENES_FILE_MODIF \
          $GT_FILE \
          $NUMBER_RESULTS \
          --force)

# Write and execute command
echo $cmd
$cmd
