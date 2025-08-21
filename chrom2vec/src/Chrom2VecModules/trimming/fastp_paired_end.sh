#!/bin/bash

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
BIND_PATH=$2

R1=$3  # Read 1 input file
R2=$4  # Read 2 input file
R1_OUT=$5  # Read 1 output file
R2_OUT=$6  # Read 2 output file
JSON_OUTPUT=$7
HTML_OUTPUT=$8

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    echo "Using singularity"
    IMAGE_URI="docker://staphb/fastp:0.23.4"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $BIND_PATH $CONTAINER  \
    fastp -i $R1 -I $R2 \
          -o $R1_OUT -O $R2_OUT \
          --thread 2 \
          -j $JSON_OUTPUT -h $HTML_OUTPUT
else
    module load fastp/0.22.0
    fastp -i $R1 -I $R2 \
          -o $R1_OUT -O $R2_OUT \
          --thread 2 \
          -j $JSON_OUTPUT -h $HTML_OUTPUT
fi
