#!/bin/bash

VIDEOS_LIST=../test_db.txt

for i in `cat $VIDEOS_LIST`; do
    in_video=../$i
    out_folder=${in_video%.*}_keyframes
    if [ ! -d "$out_folder" ]; then
        mkdir -p $out_folder
    fi
    # Get command
    cmd=$(echo python extract_keyframes.py ../$i $out_folder 1 scale=-1:480)
    # Write out command
    echo $cmd
    # Execute
    $cmd
done
