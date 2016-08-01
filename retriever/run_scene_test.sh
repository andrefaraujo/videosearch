#!/bin/bash -e

# Number of scenes to rerank
NUM_SCENES_RERANK=100

# -- Shared parameters between first and second retrieval stages
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
# Number of min words to consider a match (default: 0)
MIN_NUM_WORDS_SELECTED=0
# Word selection mode: (0=L1 norm), (1=soft assgn)
WORD_SELECTION_MODE=0

# -- Parameters for first retrieval stage
# Number of Gaussians to use in scene global descriptor
GAUSSIANS_SCENE=512
# Word selection thresh
WORD_SELECTION_THRESH_SCENE=6
# Scene parameters
SCENE_MODE=1
FPSCENE=-1

# -- Parameters for second retrieval stage (used only if NUM_SCENES_RERANK != 0)
# Number of Gaussians to use in shot global descriptor
GAUSSIANS_SHOT=512
# Word selection thresh
WORD_SELECTION_THRESH_SHOT=7
# Shot parameters
SHOT_MODE=1
FPSHOT=-1
SHOT_THRESH=0.8
GDINDEX_PATH_RERANK=../indexer/global_descriptors/test_shot_based_index_lists.sift_scfv_idx_k${GAUSSIANS_SHOT}_shot_t${SHOT_THRESH}_n${FPSHOT}_m${SHOT_MODE}
GROUP_LISTS_PATH=../indexer/global_descriptors/test_group_lists_rerank.txt

OPTIONS=""
# FIRST_FRAMES will be used for outputting correct shots/scenes
# MODE will be used for scoring correctly
if [ $NUM_SCENES_RERANK -eq 0 ]; then
    FIRST_FRAMES=../indexer/global_descriptors/test_scene_mode_${SCENE_MODE}_first_frames.txt
    MODE=$SCENE_MODE
else
    FIRST_FRAMES=../indexer/global_descriptors/test_shot_mode_${SHOT_MODE}_first_frames.txt
    MODE=$SHOT_MODE
    OPTIONS=$OPTIONS" --number_centroids_rerank "$GAUSSIANS_SHOT
    OPTIONS=$OPTIONS" --word_selection_thresh_rerank "$WORD_SELECTION_THRESH_SHOT
    OPTIONS=$OPTIONS" --gdindex_path_rerank "$GDINDEX_PATH_RERANK
    OPTIONS=$OPTIONS" --group_lists_rerank_path "$GROUP_LISTS_PATH
fi
OPTIONS=$OPTIONS" -c "$GAUSSIANS_SCENE
OPTIONS=$OPTIONS" -f "$FEAT_MODE
OPTIONS=$OPTIONS" --gdindex_parameters_path "$GDINDEX_PATH
OPTIONS=$OPTIONS" --number_output_results "$NUMBER_OUTPUT_RESULTS
OPTIONS=$OPTIONS" -v "$VERBOSE
OPTIONS=$OPTIONS" --min_number_words_visited "$MIN_NUM_WORDS_SELECTED
OPTIONS=$OPTIONS" --word_selection_mode "$WORD_SELECTION_MODE
OPTIONS=$OPTIONS" --word_selection_thresh "$WORD_SELECTION_THRESH_SCENE
OPTIONS=$OPTIONS" --shot_mode "$MODE
OPTIONS=$OPTIONS" --shot_list "$FIRST_FRAMES
OPTIONS=$OPTIONS" --number_scenes_rerank "$NUM_SCENES_RERANK

if [ $AVOID_REDUNDANT_RESULTS -eq 1 ]; then
    OPTIONS=$OPTIONS" --avoid_redundant_scene_results "
fi

# Composing output path
OUTPUT_PREFIX=test
METHOD_PATH="SCFV_scenes"
METHOD_PATH_RESULTS=${METHOD_PATH}"/gaussians_scene_"$GAUSSIANS_SCENE
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/ws_mode_"$WORD_SELECTION_MODE
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/ws_thresh_scene_"$WORD_SELECTION_THRESH_SCENE
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/min_num_words_visited_"$MIN_NUM_WORDS_SELECTED
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/scene_mode_"$SCENE_MODE
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/fpscene_"$FPSCENE
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/num_scene_rerank_"$NUM_SCENES_RERANK
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/gaussians_shot_"$GAUSSIANS_SHOT
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/ws_thresh_shot_"$WORD_SELECTION_THRESH_SHOT
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/shot_mode_"$SHOT_MODE
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/fpshot_"$FPSHOT
METHOD_PATH_RESULTS=${METHOD_PATH_RESULTS}"/shot_thresh_"$SHOT_THRESH

QUERY_LIST_FILE=test_query_list.txt
ls ../indexer/test_query/*.jpg > $QUERY_LIST_FILE
INDEX_FILE=../indexer/global_descriptors/test_scene_based_index_lists.sift_scfv_idx_k${GAUSSIANS_SCENE}_scene_n${FPSCENE}_m${SCENE_MODE}
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
