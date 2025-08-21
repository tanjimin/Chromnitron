#!/bin/bash

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
SINGULARITY_MAPPING=$2

INPUT_FASTQ_FILE=$3
OUTPUT_FASTQ_FILE=$4
JSON_OUTPUT=$5
HTML_OUTPUT=$6

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    echo "Using singularity"
    IMAGE_URI="docker://staphb/fastp:0.23.4"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $SINGULARITY_MAPPING $CONTAINER  \
    fastp -i $INPUT_FASTQ_FILE \
      -o $OUTPUT_FASTQ_FILE \
      --thread 2 \
      -j $JSON_OUTPUT -h $HTML_OUTPUT
else
    module load fastp/0.22.0
    fastp -i $INPUT_FASTQ_FILE \
      -o $OUTPUT_FASTQ_FILE \
      --thread 2 \
      -j $JSON_OUTPUT -h $HTML_OUTPUT
fi
