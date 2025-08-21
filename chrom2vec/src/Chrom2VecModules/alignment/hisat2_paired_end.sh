#!/bin/bash

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
BIND_PATH=$2
INDEX_PATH=$3
R1=$4
R2=$5
BAM_OUTPUT=$6 
SUMMARY_FILE=$7

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    IMAGE_URI="docker://biocontainers/hisat2:v2.1.0-2-deb_cv1"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    LOCAL_IMAGE_PATH=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $BIND_PATH $LOCAL_IMAGE_PATH \
    hisat2 -p 16 --mm -q --phred33 \
       -x $INDEX_PATH \
       -1 $R1 -2 $R2 \
       -S ${BAM_OUTPUT}.sam \
       --summary-file $SUMMARY_FILE

    # Sort and generate bam files
    IMAGE_URI="docker://staphb/samtools:1.17"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $BIND_PATH $CONTAINER \
    samtools view -@ 16 -bS ${BAM_OUTPUT}.sam > ${BAM_OUTPUT}.unsorted
    singularity exec --bind $BIND_PATH $CONTAINER \
    samtools sort -@ 16 -m 1G ${BAM_OUTPUT}.unsorted -o ${BAM_OUTPUT}
    singularity exec --bind $BIND_PATH $CONTAINER \
    samtools index -@ 16 ${BAM_OUTPUT}

    # Remove sam file and unsorted bam file
    rm ${BAM_OUTPUT}.sam
    rm ${BAM_OUTPUT}.unsorted

else
    module load hisat2/2.1.0
    hisat2 -p 16 --mm -q --phred33 \
       -x $INDEX_PATH \
       -1 $R1 -2 $R2 \
       -S ${BAM_OUTPUT}.sam \
       --summary-file $SUMMARY_FILE

    module load samtools/1.16
    samtools view -@ 16 -bS ${BAM_OUTPUT}.sam > ${BAM_OUTPUT}.unsorted
    samtools sort -@ 16 -m 1G ${BAM_OUTPUT}.unsorted -o ${BAM_OUTPUT}
    samtools index -@ 16 ${BAM_OUTPUT}

    rm ${BAM_OUTPUT}.sam
    rm ${BAM_OUTPUT}.unsorted
fi
