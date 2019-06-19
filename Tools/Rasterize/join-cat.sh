#!/bin/bash
#
#
#	cat style join
#
#	R. Marshall
#	July 2016
#
#	emulates MPI_Gather for output files that are partitioned row-wise
#
#
#	usage:
#	./cat-join.sh [OUTPUT_DIR [PATTERN [DELETE_OLD]]]
#
#	OUTPUT_DIR [../../output/flood2d/txt]: directory where output files reside
#	PATTERN [hOut_Mat[0-9][0-9]*_[0-9][0-9]*\.*]: regex for filename pattern
#	DELETE_OLD [0]: delete the old files, 0 to keep, nonzero int to delete
############################################################################

OUTPUT_DIR="${1:-../../output/flood2d/txt}"
WHATMAT=H
PATTERN="${2:-${WHATMAT}_[0-9][0-9]*_[0-9][0-9]*\.*}"
PATTERN0="${2:-${WHATMAT}_[0-9][0-9]*_00\.*}"
DELETE_OLD=${3:-1}

FILES=$(ls -p  $OUTPUT_DIR | grep -v / | grep -E "${PATTERN}")
PROC0_FILES=$(ls -p  $OUTPUT_DIR | grep -v / | grep -E "${PATTERN0}")
NP=0

# derive number of processes
for filename in $FILES; do
    x=$(echo $filename | cut -d'_' -f 3 | cut -d'.' -f 1)
    if [ $x -gt $NP ]; then
    	NP=$x
    fi
done

# convert MPI rank id (0-based) to number of processes (1-based)

NP=$((10#$NP + 1))

p=0

# only loop through the filenames for process 0.  We already know NP, 
#  so we can guess the remaining filenames
for filename in $PROC0_FILES; do
	BASE_NAME=$(echo $filename | cut -d'_' -f 1,2)
	echo "$BASE_NAME"
	while [ $p -lt $NP ]; do
		printf -v pp "%02d" "$p"
		cat "${OUTPUT_DIR}/${BASE_NAME}_${pp}.txt" >> "${OUTPUT_DIR}/${BASE_NAME}.txt"
		if [ $DELETE_OLD -gt 0 ]; then
			rm "${OUTPUT_DIR}/${BASE_NAME}_${pp}.txt"
		fi
		p=$((p+1))
	done
	p=0

done

exit 0
