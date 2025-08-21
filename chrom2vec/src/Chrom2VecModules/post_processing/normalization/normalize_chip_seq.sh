#!/bin/bash

MODULE_PATH=$(dirname $0)

USE_SINGULARITY=$1
BIND_PATH=$2

PEAK_FILE=$3
NON_PEAK_FILE=$4
CHIP_IP_ZARR=$5
CHIP_INPUT_ZARR=$6
OUTPUT_ZARR=$7

# Check if both PEAK_FILE and NON_PEAK_FILE exist
if [ ! -f "$PEAK_FILE" ] || [ ! -f "$NON_PEAK_FILE" ]; then
    echo "One or both of the peak files do not exist. Downloading..."
    bash $MODULE_PATH/download_atac_seq_peaks.sh $USE_SINGULARITY $BIND_PATH $PEAK_FILE $NON_PEAK_FILE
fi

# If use singularity
if [ $USE_SINGULARITY == "True" ]; then
    echo "Using singularity"
    IMAGE_URI="docker://tanjimin/pybigwig:0.0.1"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $BIND_PATH $CONTAINER \
    python3 $MODULE_PATH/normalize_chip_seq.py $PEAK_FILE $NON_PEAK_FILE $CHIP_IP_ZARR $CHIP_INPUT_ZARR $OUTPUT_ZARR
else
    source ~/.bashrc
    conda activate /gpfs/data/tsirigoslab/home/jt3545/conda/hic
    python3 $MODULE_PATH/normalize_chip_seq.py $PEAK_FILE $NON_PEAK_FILE $CHIP_IP_ZARR $CHIP_INPUT_ZARR $OUTPUT_ZARR
fi

