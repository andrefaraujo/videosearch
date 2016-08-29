#!/bin/bash -e

# This script extracts HesAff keypoints + SIFT descriptors, as
# described in http://lear.inrialpes.fr/~jegou/data.php
# The $COMPUTE_DESCRIPTORS_CMD program can be downloaded on the
# same website; make sure it is located in this directory

# Check platform in use
if [ "$(uname)" == "Darwin" ]; then
    # OS X
    COMPUTE_DESCRIPTORS_CMD=compute_descriptors_mac
else
    # Linux
    COMPUTE_DESCRIPTORS_CMD=compute_descriptors_linux64
fi

# Get list of keyframes
IMAGE_LIST=test_all_db_keyframes.txt
ls ../test_db/*/*.jpg > $IMAGE_LIST

# Change this path to where netpbm binaries are located
NETPBM_BIN_PATH=/usr/local/bin/

for i in `cat $IMAGE_LIST`; do

    # Renaming
    infile=$i
    tmpfile=${infile/jpg/pgm}
    outfile=${infile/jpg/siftgeo}

    # Check if output file already exists, if so we skip it
	if [ -f "$outfile" ]; then
		continue
	fi

    # Rescaling and intensity normalization 
    cmd_pgm_1=$(echo djpeg $infile)
    cmd_pgm_2=$(echo $NETPBM_BIN_PATH/ppmtopgm)
    cmd_pgm_3=$(echo $NETPBM_BIN_PATH/pnmnorm -bpercent=0.01 -wpercent=0.01 -maxexpand=400)
    cmd_pgm_4=$(echo $NETPBM_BIN_PATH/pamscale -pixels $[1024*768])
    
    # Write and execute command
    echo "time $cmd_pgm_1 | $cmd_pgm_2 | $cmd_pgm_3 | $cmd_pgm_4 > $tmpfile"
    time $cmd_pgm_1 | $cmd_pgm_2 | $cmd_pgm_3 | $cmd_pgm_4 > $tmpfile

    # Compute descriptors
    cmd_compute=$(echo time \
        ./${COMPUTE_DESCRIPTORS_CMD} \
        -i $tmpfile \
        -o4 $outfile \
        -hesaff \
        -sift)

    # Write and execute command
    echo $cmd_compute
    $cmd_compute
done
