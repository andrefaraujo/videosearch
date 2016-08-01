#!/bin/bash -e

# Get list of queries
IMAGE_LIST=test_query_list.txt
ls ../test_query/*.jpg > $IMAGE_LIST

# Choose the number of threads
NUMBER_THREADS=1

./extract -t ${NUMBER_THREADS} -i ${IMAGE_LIST} -s
