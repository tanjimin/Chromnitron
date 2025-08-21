#!/bin/bash

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
SINGULARITY_MAPPING=$2

INDEX_PATH=$3
READ_FILE=$4
OUTPUT_BAM_FILE=$5
SUMMARY_FILE=$6

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    IMAGE_URI="docker://biocontainers/hisat2:v2.1.0-2-deb_cv1"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    LOCAL_IMAGE_PATH=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $SINGULARITY_MAPPING $LOCAL_IMAGE_PATH \
    hisat2 -p 16 --mm -q --phred33 \
       -x $INDEX_PATH \
       -U $READ_FILE \
       -S $OUTPUT_BAM_FILE.sam \
       --summary-file $SUMMARY_FILE

    # Sort and generate bam files
    IMAGE_URI="docker://staphb/samtools:1.17"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $SINGULARITY_MAPPING $CONTAINER \
    samtools view -@ 16 -bS $OUTPUT_BAM_FILE.sam > $OUTPUT_BAM_FILE.unsorted
    singularity exec --bind $SINGULARITY_MAPPING $CONTAINER \
    samtools sort -@ 16 -m 1G $OUTPUT_BAM_FILE.unsorted -o $OUTPUT_BAM_FILE

    # Index bam files
    singularity exec --bind $SINGULARITY_MAPPING $CONTAINER \
    samtools index -@ 16 $OUTPUT_BAM_FILE

    # Remove sam file and unsorted bam file
    rm $OUTPUT_BAM_FILE.sam
    rm $OUTPUT_BAM_FILE.unsorted

else
    module load hisat2/2.1.0
    hisat2 -p 16 --mm -q --phred33 \
       -x $INDEX_PATH \
       -U $READ_FILE \
       -S $OUTPUT_BAM_FILE.sam \
       --summary-file $SUMMARY_FILE

    module load samtools/1.16
    samtools view -@ 16 -bS $OUTPUT_BAM_FILE.sam > $OUTPUT_BAM_FILE.unsorted
    samtools sort -@ 16 -m 1G $OUTPUT_BAM_FILE.unsorted -o $OUTPUT_BAM_FILE
    samtools index -@ 16 $OUTPUT_BAM_FILE

    rm $OUTPUT_BAM_FILE.sam
    rm $OUTPUT_BAM_FILE.unsorted
fi
