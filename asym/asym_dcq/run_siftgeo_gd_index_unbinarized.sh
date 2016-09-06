#!/bin/bash -e

# Parameters for index
GDINDEX_PATH=../../indexer/global_descriptors/trained_parameters
CENTROIDS=512
LD_MODE=1
VERBOSE=1

# Database file 
list=database_list_public.txt

if [ $LD_MODE -eq 0 ]; then
    LD_NAME=sift
elif [ $LD_MODE -eq 1 ]; then
    LD_NAME=siftgeo
else
    echo "Unrecognized LD_MODE"
    exit
fi

# Compose output index name
out_index=${list%.txt}.${LD_NAME}_fv_idx_k${CENTROIDS}

# Command line
cmd=$(echo time \
    ../../indexer/global_descriptors/index_dataset \
    -i $list \
    -o $out_index \
    -r $GDINDEX_PATH \
    -c $CENTROIDS \
    -l $LD_MODE \
    -v $VERBOSE \
    --gd_unbinarized \
    --gd_intra_normalization)

# Write and execute command
echo $cmd
$cmd

