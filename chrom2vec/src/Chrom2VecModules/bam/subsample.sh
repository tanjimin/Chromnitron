#!/bin/bash

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
SINGULARITY_MAPPING=$2

RANDOM_SEED=$3
SUBSAMPLE_READS=$4
INPUT_BAM=$5
OUTPUT_BAM=$6

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    IMAGE_URI="docker://staphb/samtools:1.17"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif

    # Total number of reads
    TOTAL_READS=$(singularity exec --bind $SINGULARITY_MAPPING $CONTAINER \
        samtools view -@ 16 -c $INPUT_BAM)
    echo "Total number of reads: $TOTAL_READS"

    # Calculate subsample fraction
    FRACTION=$(echo "scale=20; ${SUBSAMPLE_READS} / ${TOTAL_READS}" | bc)
    echo "Fraction of reads to keep: $FRACTION"

    # Check if fraction is greater than 1
    if (( $(echo "$FRACTION > 1" | bc -l) )); then
        FRACTION=1
        echo "WARNING: Fraction of reads to keep is greater than 1. CONTINUE."
        rm -f $OUTPUT_BAM
        cp $INPUT_BAM $OUTPUT_BAM
        rm -f $OUTPUT_BAM.bai
        cp $INPUT_BAM.bai $OUTPUT_BAM.bai
        exit 0
    fi

    # Subsample bam file
    singularity exec --bind $SINGULARITY_MAPPING $CONTAINER \
    samtools view -@ 16 -bs $RANDOM_SEED$FRACTION -o $OUTPUT_BAM $INPUT_BAM

    # Index bam files
    singularity exec --bind $SINGULARITY_MAPPING $CONTAINER \
    samtools index -@ 16 $OUTPUT_BAM

else
    module load samtools/1.16
    TOTAL_READS=$(samtools view -@ 16 -c $INPUT_BAM)
    echo "Total number of reads: $TOTAL_READS"
    FRACTION=$(echo "scale=20; ${SUBSAMPLE_READS} / ${TOTAL_READS}" | bc)
    echo "Fraction of reads to keep: $FRACTION"
    if (( $(echo "$FRACTION > 1" | bc -l) )); then
        FRACTION=1
        echo "WARNING: Fraction of reads to keep is greater than 1. CONTINUE."
        rm -f $OUTPUT_BAM
        cp $INPUT_BAM $OUTPUT_BAM
        rm -f $OUTPUT_BAM.bai
        cp $INPUT_BAM.bai $OUTPUT_BAM.bai
        exit 0
    fi
    samtools view -@ 16 -bs $RANDOM_SEED$FRACTION -o $OUTPUT_BAM $INPUT_BAM
    samtools index -@ 16 $OUTPUT_BAM
fi

