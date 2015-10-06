#!/bin/bash

# Number of Gaussians to use in global descriptor
GAUSSIANS=512
# Feature mode: 0=SIFT; others currently not supported
FEAT_MODE=0
# Path where to find GDIndex's trained parameters
GDINDEX_PATH=../indexer/global_descriptors/trained_parameters
# Number of results to output
NUMBER_OUTPUT_RESULTS=100
# Flag that if set avoids outputting redundant scene results
AVOID_REDUNDANT_RESULTS=1
# Verbose level
VERBOSE=1
# Number of min words to consider a match (default: 20)
MIN_NUM_WORDS_SELECTED=0
# Word selection mode: (0=L1 norm), (1=soft assgn)
WORD_SELECTION_MODE=0
# Word selection thresh
WORD_SELECTION_THRESH=8

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
OUTPUT_PREFIX=test
METHOD_PATH="SCFV_frames"
METHOD_PATH_RESULTS=${METHOD_PATH}"/gaussians_"$GAUSSIANS
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/ws_mode_"$WORD_SELECTION_MODE
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/ws_thresh_"$WORD_SELECTION_THRESH
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/min_num_words_visited_"$MIN_NUM_WORDS_SELECTED

QUERY_LIST_FILE=test_query_list.txt
ls ../indexer/test_query/*.jpg > $QUERY_LIST_FILE
INDEX_FILE=../indexer/global_descriptors/test_frame_based_index_lists.sift_scfv_idx_k512
DB_LIST_FILE=test_frame_db_list.txt
ls ../indexer/test_db/*/*.jpg > $DB_LIST_FILE
OUTPUT_PATH=$OUTPUT_PREFIX/results/${METHOD_PATH_RESULTS}
OUTPUT_BASE=${OUTPUT_PATH}/out

# Create output folder if necessary
if [ ! -d "$OUTPUT_PATH" ]; then
	mkdir -p $OUTPUT_PATH
fi

# Compose command line
cmd=$(echo \
	./retrieve_on_dataset \
	-i $INDEX_FILE \
	-d $DB_LIST_FILE \
	-q $QUERY_LIST_FILE \
	-o $OUTPUT_BASE \
	$OPTIONS)

# Write and execute command
echo $cmd
$cmd
