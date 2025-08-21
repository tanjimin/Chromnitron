#!/bin/bash
#
# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
SINGULARITY_MAPPING=$2

OUTPUT_BAM=$3
# Assuming the rest of the input files start from $4 onwards
INPUT_FILES=${@:4}

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    IMAGE_URI="docker://staphb/samtools:1.17"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    # Merge all bam files
    singularity exec --bind $SINGULARITY_MAPPING $CONTAINER \
    samtools merge -@ 16 -o $OUTPUT_BAM $INPUT_FILES
    # Index bam files
    singularity exec --bind $SINGULARITY_MAPPING $CONTAINER \
    samtools index -@ 16 $OUTPUT_BAM
else
    module load samtools/1.16
    samtools merge -@ 16 -o $OUTPUT_BAM $INPUT_FILES
    samtools index -@ 16 $OUTPUT_BAM
fi

