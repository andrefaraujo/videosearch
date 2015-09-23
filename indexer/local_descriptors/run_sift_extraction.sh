#!/bin/bash

IMAGE_LIST=test_images.txt

# Choose the number of threads
NUMBER_THREADS=5

./extract -t ${NUMBER_THREADS} -i ${IMAGE_LIST} -s
