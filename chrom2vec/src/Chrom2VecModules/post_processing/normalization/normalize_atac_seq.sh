#!/bin/bash

MODULE_PATH=$(dirname $0)

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
BIND_PATH=$2

PEAK_FILE=$3
NON_PEAK_FILE=$4
ZARR_TRACK_FILE=$5
OUTPUT_FILE=$6

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    echo "Using singularity"
    IMAGE_URI="docker://tanjimin/pybigwig:0.0.1"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $BIND_PATH $CONTAINER \
    python3 $MODULE_PATH/normalize_atac_seq.py $PEAK_FILE $NON_PEAK_FILE $ZARR_TRACK_FILE $OUTPUT_FILE
else
    source ~/.bashrc
    conda activate /gpfs/data/tsirigoslab/home/jt3545/conda/hic
    python3 $MODULE_PATH/normalize_atac_seq.py $PEAK_FILE $NON_PEAK_FILE $ZARR_TRACK_FILE $OUTPUT_FILE
fi
