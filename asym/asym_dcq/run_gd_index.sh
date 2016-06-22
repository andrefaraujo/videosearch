#!/bin/bash

# Parameters for index
GDINDEX_PATH=../../indexer/global_descriptors/trained_parameters
CENTROIDS=512
LD_MODE=0
VERBOSE=1

# Database file 
list=database_list_public.txt

# Compose output index name
out_index=${list%.txt}.sift_scfv_idx_k${CENTROIDS}

# Command line
cmd=$(echo time \
    ../../indexer/global_descriptors/index_dataset \
    -i $list \
    -o $out_index \
    -r $GDINDEX_PATH \
    -c $CENTROIDS \
    -l $LD_MODE \
    -v $VERBOSE)

# Write and execute command
echo $cmd
$cmd

