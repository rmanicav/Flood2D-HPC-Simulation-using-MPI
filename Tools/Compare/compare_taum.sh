#!/bin/bash

# ./compare_taum.sh 1136 624 60 ../../output/flood2d-output-pkg-1/txt ../../output/flood2d-output-pkg-2/txt ../../output/flood2d-compare

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)


NROWS="$1"    # num rows
NCOLS="$2"    # num cols
NFILES="$3"   # num files
DIR_A="$4"    # dir of input files
DIR_B="$5"    # dir of input files
DIR_C="$6"    # dir of output files


BGFILE="${PROJECT_ROOT}/input/dem/ts_dam_break.dem"


rm -rf "$DIR_C"
mkdir -p "$DIR_C"

for n in $(seq 1 "$NFILES"); do
	printf -v fnum "%02d" "$n"
	./compare "$NROWS" "$NCOLS" "${DIR_A}/H_${fnum}.txt" "${DIR_B}/H_${fnum}.txt" >"${DIR_C}/H_${fnum}.txt"

done
