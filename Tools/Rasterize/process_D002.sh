#!/bin/bash

START=$PWD
PROJECT_ROOT=$(cd ../../ && echo $PWD && cd $START)
OUTPUT_ROOT="${1:-../../output}"
OUTPUT_NAME="${2:-flood2d-output-pkg}"
PREFIX="${3:-H}"
CFGFILE="${4:-${PROJECT_ROOT}/input/cfg/D002.cfg}"
HYDROFILE="${5:-${PROJECT_ROOT}/input/hyg/D002.txt}"
DESC="${6:-D002 Simulation}"
INPUT_DIR="${PROJECT_ROOT}/output/flood2d/txt"
BGFILE="${PROJECT_ROOT}/input/dem/D002.txt"

#set -x

rm -rf "$OUTPUT_ROOT/$OUTPUT_NAME" #2>/dev/null
mkdir "$OUTPUT_ROOT/$OUTPUT_NAME"
mkdir "$OUTPUT_ROOT/$OUTPUT_NAME/txt"
mkdir "$OUTPUT_ROOT/$OUTPUT_NAME/gif"


./join-cat.sh "$INPUT_DIR"
mv ${INPUT_DIR}/${PREFIX}_*.txt ${OUTPUT_ROOT}/${OUTPUT_NAME}

echo "calling rasterize"

./rasterize --datadir="${PROJECT_ROOT}/output/${OUTPUT_NAME}" --nrows=2394 --ncols=3159 --bgcolors=128 --fgcolors=128 --bgdata="${BGFILE}"
 
echo "moving to package folder"
mv ${PROJECT_ROOT}/output/${OUTPUT_NAME}/${PREFIX}_*.gif ${PROJECT_ROOT}/output/${OUTPUT_NAME}/gif
mv ${PROJECT_ROOT}/output/${OUTPUT_NAME}/${PREFIX}_*.txt ${PROJECT_ROOT}/output/${OUTPUT_NAME}/txt

cp "${CFGFILE}" "${PROJECT_ROOT}/output/${OUTPUT_NAME}"
cp "${HYDROFILE}" "${PROJECT_ROOT}/output/${OUTPUT_NAME}"

rm -rf "${PROJECT_ROOT}/output/${OUTPUT_NAME}/"*.bmp

echo "$DESC"
echo $(md5sum "${BGFILE}") "${BGFILE}"
echo "${INPUT_DIR}"
echo "${BGFILE}"

cd $START

