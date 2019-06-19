#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
INPUT_DIR="${PROJECT_ROOT}/output/flood2d/bin/"
OUTPUT_DIR="${PROJECT_ROOT}/output/flood2d/txt/"

./bin2ascii $INPUT_DIR $OUTPUT_DIR 7183 9478 1

cd $START