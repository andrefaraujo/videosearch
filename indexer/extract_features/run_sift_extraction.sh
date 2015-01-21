#!/bin/bash

# Note: sometimes there is an error when running this, such as "./extract: error while loading shared libraries: libvl.so: cannot open shared object file: No such file or directory"
# To resolve it, do something like:
# sudo cp ../../common/vlfeat-0.9.18/bin/glnxa64/libvl.so /usr/lib/

IMAGE_LIST=test_images.txt

# Choose the number of threads
NUMBER_THREADS=5

./extract -t ${NUMBER_THREADS} -i ${IMAGE_LIST} -s
