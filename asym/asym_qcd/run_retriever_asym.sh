#!/bin/bash -e

# Parameter that defines the clutter level
C=5
C_printf=$(printf "%04d" $C)
# Number of Gaussians to use in global descriptor
GAUSSIANS=512
# Feature mode: 0=SIFT; others currently not supported
FEAT_MODE=0
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
WORD_SELECTION_THRESH=0.0004
# "Shot" parameters
SHOT_MODE=1
FPSHOT=-1

OPTIONS=""
OPTIONS=$OPTIONS" -c "$GAUSSIANS
OPTIONS=$OPTIONS" -f "$FEAT_MODE
OPTIONS=$OPTIONS" --gdindex_parameters_path "$GDINDEX_PATH
OPTIONS=$OPTIONS" --number_output_results "$NUMBER_OUTPUT_RESULTS
OPTIONS=$OPTIONS" -v "$VERBOSE
OPTIONS=$OPTIONS" --min_number_words_visited "$MIN_NUM_WORDS_SELECTED
OPTIONS=$OPTIONS" --word_selection_mode "$WORD_SELECTION_MODE
OPTIONS=$OPTIONS" --word_selection_thresh "$WORD_SELECTION_THRESH

if [ $AVOID_REDUNDANT_RESULTS -eq 1 ]; then
    OPTIONS=$OPTIONS" --avoid_redundant_scene_results "
fi

# Composing output path
OUTPUT_PREFIX=results/"C_"${C_printf}
METHOD_PATH="SCFV_shots"
METHOD_PATH_RESULTS=${METHOD_PATH}"/gaussians_"$GAUSSIANS
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/ws_mode_"$WORD_SELECTION_MODE
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/ws_thresh_"$WORD_SELECTION_THRESH
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/min_num_words_visited_"$MIN_NUM_WORDS_SELECTED

QUERY_LIST_FILE=query_list_public.txt
INDEX_FILE=db_list_clutter${C_printf}_public.sift_scfv_idx_k${GAUSSIANS}_shot_n${FPSHOT}_m${SHOT_MODE}
DB_LIST_FILE=db_list_clutter0000_public.txt # list with uncluttered DB images
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
	-q $QUERY_LIST_FILE \
	-o $OUTPUT_BASE \
	$OPTIONS)

# Write and execute command
echo $cmd
$cmd
