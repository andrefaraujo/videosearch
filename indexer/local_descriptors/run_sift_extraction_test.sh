#!/bin/bash -e

# Get list of keyframes
IMAGE_LIST=test_all_db_keyframes.txt
ls ../test_db/*/*.jpg > $IMAGE_LIST

# Choose the number of threads
NUMBER_THREADS=5

./extract -t ${NUMBER_THREADS} -i ${IMAGE_LIST} -s
