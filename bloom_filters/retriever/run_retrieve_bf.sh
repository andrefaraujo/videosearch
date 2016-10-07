#!/bin/bash -e

# ----- Mandatory arguments
# Path to list of video clips
CLIP_LIST_PATH=clip_list_path.txt
ls ../../indexer/test_db/*.mp4 > $CLIP_LIST_PATH
# Path to file with list of queries
QUERY_LIST_PATH=../../indexer/global_descriptors/test_query_list.txt
# Number of bits per hash
NUM_BITS=16
# Number of hashers
NUM_HASHERS=512

# ----- Options
# Number of Gaussians to use in point-indexed
GAUSSIANS=512
# Feature mode: 0=SIFT; 1=SIFTGeo
FEAT_MODE=1
# Number of results to output
RESULTS_PER_QUERY=100
# TF-IDF power normalization
ALPHA=0.75
# Verbose level
VERBOSE=4

if [ $FEAT_MODE -eq 0 ]; then
    FEAT_NAME=sift
elif [ $FEAT_MODE -eq 1 ]; then
    FEAT_NAME=siftgeo
else
    echo "Unrecognized FEAT_MODE"
    exit
fi

OPTIONS=""
OPTIONS=$OPTIONS" --number_gaussians "$GAUSSIANS
OPTIONS=$OPTIONS" --feature_mode "$FEAT_MODE
OPTIONS=$OPTIONS" --results_per_query "$RESULTS_PER_QUERY
OPTIONS=$OPTIONS" --alpha "$ALPHA
OPTIONS=$OPTIONS" --verbose_level "$VERBOSE

# Composing output path
OUTPUT_PREFIX=test_bf_results
METHOD_PATH_RESULTS=${METHOD_PATH}"/feat_name_"$FEAT_NAME
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/gaussians_"$GAUSSIANS
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/alpha_"$ALPHA
OUTPUT_PATH=$OUTPUT_PREFIX/results/${METHOD_PATH_RESULTS}

# Create output folder if necessary
if [ ! -d "$OUTPUT_PATH" ]; then
	mkdir -p $OUTPUT_PATH
fi

# Compose command line
cmd=$(echo \
	      ./retrieve_on_dataset_bf \
	      --list_path $CLIP_LIST_PATH \
          --query_list $QUERY_LIST_PATH \
          --output $OUTPUT_PATH \
          --n_bits $NUM_BITS \
          --n_hashers $NUM_HASHERS \
	      $OPTIONS)

# Write and execute command
echo $cmd
$cmd
