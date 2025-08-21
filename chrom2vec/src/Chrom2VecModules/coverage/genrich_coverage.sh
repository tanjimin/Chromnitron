#!/bin/bash

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
BIND_PATH=$2

INPUT_BAM_FILE=$3 
BLACKLIST_REGIONS=$4 
ASSEMBLY_NAME=$5
OUTPUT_BEDGRAPH_FILE=$6

NAME_SORTED_BAM_FILE=${INPUT_BAM_FILE}.name_sorted.bam

# Check if BLACKLIST_REGIONS exists, download if it doesn't
if [ ! -f "$BLACKLIST_REGIONS" ]; then
    echo "Blacklist file not found. Downloading..."
    wget https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/${ASSEMBLY_NAME}-blacklist.v2.bed.gz -O ${BLACKLIST_REGIONS}.gz
    gunzip ${BLACKLIST_REGIONS}.gz
fi

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    IMAGE_URI="docker://staphb/samtools:1.17"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $BIND_PATH $CONTAINER \
    samtools sort -n -@ 8 -m 1G $INPUT_BAM_FILE -o $NAME_SORTED_BAM_FILE

    IMAGE_URI="docker://juettemann/genrich:0.6"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $BIND_PATH $CONTAINER \
    Genrich -j -y -r -X \
        -t $NAME_SORTED_BAM_FILE \
        -e chrM \
        -E $BLACKLIST_REGIONS \
        -d 20 \
        -m 30 \
        -k ${OUTPUT_BEDGRAPH_FILE}.with_header

    singularity exec --bind $BIND_PATH $CONTAINER \
    tail -n +3 ${OUTPUT_BEDGRAPH_FILE}.with_header | awk '{ print $1"\t"$2"\t"$3"\t"$4 }' > $OUTPUT_BEDGRAPH_FILE
else
    module load samtools/1.16
    samtools sort -n -@ 8 -m 1G $INPUT_BAM_FILE -o $NAME_SORTED_BAM_FILE
    module load genrich/0.6
    Genrich -j -y -r -X \
        -t $INPUT_BAM_FILE \
        -e chrM \
        -E $BLACKLIST_REGIONS \
        -d 20 \
        -m 30 \
        -k ${OUTPUT_BEDGRAPH_FILE}.with_header \ &&
    tail -n +3 ${OUTPUT_BEDGRAPH_FILE}.with_header | awk '{ print $1"\t"$2"\t"$3"\t"$4 }' > $OUTPUT_BEDGRAPH_FILE
fi
