#!/bin/bash

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
BIND_PATH=$2

BAM_FILE=$3
BLACKLIST_FILE=$4
ASSEMBLY_NAME=$5
OUTPUT_FILE=$6

# Check if BLACKLIST_FILE exists, download if it doesn't
if [ ! -f "$BLACKLIST_FILE" ]; then
    echo "Blacklist file not found. Downloading..."
    wget https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/${ASSEMBLY_NAME}-blacklist.v2.bed.gz -O ${BLACKLIST_FILE}.gz
    gunzip ${BLACKLIST_FILE}.gz
fi

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    IMAGE_URI="docker://stjudecloud/deeptools:branch-encoding-1.0.0"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $BIND_PATH $CONTAINER \
    bamCoverage -b $BAM_FILE \
                -of bedgraph \
                --blackListFileName $BLACKLIST_FILE \
                -o $OUTPUT_FILE \
                --binSize 1 \
                --ignoreDuplicates \
                --minMappingQuality 30 \
                -p 16
else
    module load deeptools/3.5.1
    bamCoverage -b $BAM_FILE \
                -of bedgraph \
                --blackListFileName $BLACKLIST_FILE \
                -o $OUTPUT_FILE \
                --binSize 1 \
                --ignoreDuplicates \
                --minMappingQuality 30 \
                -p 16
fi

