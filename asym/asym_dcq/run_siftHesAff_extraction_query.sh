#!/bin/bash -e

# Check platform in use
if [ "$(uname)" == "Darwin" ]; then
    # OS X
    COMPUTE_DESCRIPTORS_CMD=compute_descriptors_mac
else
    # Linux
    COMPUTE_DESCRIPTORS_CMD=compute_descriptors_linux64
fi

CMD_PATH=../../indexer/local_descriptors/

# Get list of images
IMAGE_LIST=query_list_clutter0005_public.txt

# Change this path to where netpbm binaries are located
NETPBM_BIN_PATH=/usr/local/netpbm/bin/

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
        ${CMD_PATH}/${COMPUTE_DESCRIPTORS_CMD} \
        -i $tmpfile \
        -o4 $outfile \
        -hesaff \
        -sift)

    # Write and execute command
    echo $cmd_compute
    $cmd_compute
done
