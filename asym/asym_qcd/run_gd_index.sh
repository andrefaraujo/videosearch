#!/bin/bash

# Parameter that defines the clutter level
C=5
C_printf=$(printf "%04d" $C)

# Parameters for index
GDINDEX_PATH=../../indexer/global_descriptors/trained_parameters
CENTROIDS=512
LD_MODE=0
VERBOSE=1

# "Shot" parameters -- used to extract FV from a set of images
SHOT_MODE=1
SHOT_KEYF=-1

# Database file 
list=db_list_clutter${C_printf}_public.txt

# Inds file (defines which sets of images are grouped together)
shot_res_file=inds_list_clutter${C_printf}_public.txt

# Compose output index name
out_index=${list%.txt}.sift_scfv_idx_k${CENTROIDS}_shot_n${SHOT_KEYF}_m${SHOT_MODE}

# Command line
cmd=$(echo time \
    ../../indexer/global_descriptors/index_dataset \
    -i $list \
    -o $out_index \
    -r $GDINDEX_PATH \
    -c $CENTROIDS \
    -l $LD_MODE \
    -v $VERBOSE \
    -m $SHOT_MODE \
    -k $SHOT_KEYF \
    -s $shot_res_file)

# Write and execute command
echo $cmd
$cmd

