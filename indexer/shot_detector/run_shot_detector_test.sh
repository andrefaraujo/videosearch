#!/bin/bash

THRESHOLD=0.8
VERBOSE=0

# Get list of clips (make sure to extract keyframes before running this)
lists=test_shot_lists.txt
ls ../test_db/*.txt > $lists

for list in `cat $lists`
do
    out_name=${list%.txt}.shot_t${THRESHOLD}

    # Get command
    cmd=$(echo ./shot_detector -l $list -t $THRESHOLD \
        -o $out_name -v $VERBOSE)
    # Write out command
    echo $cmd
    # Run command
    $cmd
done