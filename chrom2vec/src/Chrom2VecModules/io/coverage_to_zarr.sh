#!/bin/bash

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
BIND_PATH=$2
INPUT_BEDGRAPH=$3  # Replace INPUT_BEDGRAPH with a more descriptive name based on its purpose
OUTPUT_ZARR=$4  # Replace OUTPUT_ZARR with a more descriptive name based on its purpose

# Determine the module path
MODULE_PATH=$(dirname $0)

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    echo "Using singularity"
    IMAGE_URI="docker://tanjimin/seq2vec:0.0.1"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $BIND_PATH $CONTAINER \
    python3 $MODULE_PATH/coverage_to_zarr.py $INPUT_BEDGRAPH $OUTPUT_ZARR
else
    source ~/.bashrc
    conda activate /gpfs/data/tsirigoslab/home/jt3545/conda/corigami
    python3 $MODULE_PATH/coverage_to_zarr.py $INPUT_BEDGRAPH $OUTPUT_ZARR
fi
