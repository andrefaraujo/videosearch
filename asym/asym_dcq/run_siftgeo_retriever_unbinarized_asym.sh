#!/bin/bash -e

# Parameter that defines the clutter level
C=5
C_printf=$(printf "%04d" $C)
# Number of Gaussians to use in global descriptor
GAUSSIANS=512
# Feature mode: 0=SIFT; 1=SIFTGEO; others currently not supported
FEAT_MODE=1
# Path where to find GDIndex's trained parameters
GDINDEX_PATH=../../indexer/global_descriptors/trained_parameters
# Number of results to output
NUMBER_OUTPUT_RESULTS=100
# Flag that if set avoids outputting redundant scene results
AVOID_REDUNDANT_RESULTS=0
# Verbose level
VERBOSE=1
# Number of min words to consider a match (default: 0)
MIN_NUM_WORDS_SELECTED=0
# Word selection mode: (0=L1 norm), (1=soft assgn)
WORD_SELECTION_MODE=1
# Word selection thresh
WORD_SELECTION_THRESH=0.0016
# "Shot" parameters
SHOT_MODE=1
FPSHOT=-1
# Asym. scoring mode
ASYM_SCORING_MODE=2
# Score denominator power
SCORE_DEN_POWER_NORM=-0.1

if [ $FEAT_MODE -eq 0 ]; then
    FEAT_NAME=sift
elif [ $FEAT_MODE -eq 1 ]; then
    FEAT_NAME=siftgeo
else
    echo "Unrecognized FEAT_MODE"
    exit
fi

OPTIONS=""
OPTIONS=$OPTIONS" --gd_unbinarized"
OPTIONS=$OPTIONS" --gd_intra_normalization"
OPTIONS=$OPTIONS" -c "$GAUSSIANS
OPTIONS=$OPTIONS" -f "$FEAT_MODE
OPTIONS=$OPTIONS" --gdindex_parameters_path "$GDINDEX_PATH
OPTIONS=$OPTIONS" --number_output_results "$NUMBER_OUTPUT_RESULTS
OPTIONS=$OPTIONS" -v "$VERBOSE
OPTIONS=$OPTIONS" --min_number_words_visited "$MIN_NUM_WORDS_SELECTED
OPTIONS=$OPTIONS" --word_selection_mode "$WORD_SELECTION_MODE
OPTIONS=$OPTIONS" --word_selection_thresh "$WORD_SELECTION_THRESH
OPTIONS=$OPTIONS" --score_den_power_norm "$SCORE_DEN_POWER_NORM
OPTIONS=$OPTIONS" --asym_scoring_mode "$ASYM_SCORING_MODE

if [ $AVOID_REDUNDANT_RESULTS -eq 1 ]; then
    OPTIONS=$OPTIONS" --avoid_redundant_scene_results "
fi

# Composing output path
OUTPUT_PREFIX=results/"C_"${C_printf}
METHOD_PATH="FV_images"
METHOD_PATH_RESULTS=${METHOD_PATH}"/gaussians_"$GAUSSIANS
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/ws_mode_"$WORD_SELECTION_MODE
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/ws_thresh_"$WORD_SELECTION_THRESH
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/score_den_power_"$SCORE_DEN_POWER_NORM
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/min_num_words_visited_"$MIN_NUM_WORDS_SELECTED
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/asym_scoring_mode_"$ASYM_SCORING_MODE

QUERY_LIST_FILE=query_list_clutter0000_public.txt # list with uncluttered query images
QUERY_CLUTTERED_FILE=query_list_clutter${C_printf}_public.txt
QUERY_INDEX_FILE=${QUERY_CLUTTERED_FILE%.txt}.${FEAT_NAME}_fv_idx_k${GAUSSIANS}_shot_n${FPSHOT}_m${SHOT_MODE}
INDEX_FILE=database_list_public.${FEAT_NAME}_fv_idx_k${GAUSSIANS}
DB_LIST_FILE=database_list_public.txt
OUTPUT_PATH=$OUTPUT_PREFIX/${METHOD_PATH_RESULTS}
OUTPUT_BASE=${OUTPUT_PATH}/out

# Create output folder if necessary
if [ ! -d "$OUTPUT_PATH" ]; then
    mkdir -p $OUTPUT_PATH
fi

# Compose command line
cmd=$(echo \
    ../../retriever/retrieve_on_dataset \
    -i $INDEX_FILE \
    -d $DB_LIST_FILE \
    --query_index $QUERY_INDEX_FILE \
    -q $QUERY_LIST_FILE \
    -o $OUTPUT_BASE \
    $OPTIONS)

# Write and execute command
echo $cmd
$cmd
