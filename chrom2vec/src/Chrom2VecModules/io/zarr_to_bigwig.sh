#!/bin/bash

MODULE_PATH=$(dirname $0)

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
BIND_PATH=$2

INPUT_ZARR=$3  # Replace INPUT_ZARR with a more descriptive name based on its purpose
OUTPUT_BIGWIG=$4  # Replace OUTPUT_BIGWIG with a more descriptive name based on its purpose

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    echo "Using singularity"
    IMAGE_URI="docker://tanjimin/pybigwig:0.0.1"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $BIND_PATH $CONTAINER \
    python3 $MODULE_PATH/zarr_to_bigwig.py $INPUT_ZARR $OUTPUT_BIGWIG
else
    source ~/.bashrc
    conda activate /gpfs/data/tsirigoslab/home/jt3545/conda/hic
    python3 $MODULE_PATH/zarr_to_bigwig.py $INPUT_ZARR $OUTPUT_BIGWIG
fi

