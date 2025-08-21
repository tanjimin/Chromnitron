#!/bin/bash

USE_SINGULARITY=$1
SINGULARITY_MAPPING=$2

INPUT_FATSTQ_FILE=$3
OUTPUT_DIR=$4

mkdir -p $OUTPUT_DIR

if [ $USE_SINGULARITY == "True" ]; then
    echo "Using singularity"
    IMAGE_URI="docker://staphb/fastqc:0.12.1"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $SINGULARITY_MAPPING $CONTAINER  \
    fastqc $INPUT_FATSTQ_FILE -o $OUTPUT_DIR \
        -t 16
else
    module load fastqc/0.11.7
    fastqc $INPUT_FATSTQ_FILE -o $OUTPUT_DIR \
        -t 16
fi
