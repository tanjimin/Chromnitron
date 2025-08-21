#!/bin/bash

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
BIND_PATH=$2

ASSEMBLY=$3
PEAK_TYPE=$4
IP_SAMPLE=$5
INPUT_SAMPLE=$6
CHR_SIZES=$7
OUTPUT_PEAK=$8
OUTPUT_NON_PEAK=$9

OUTPUT_PATH=$(dirname $OUTPUT_PEAK)

# Check if CHR_SIZES file exists, download if it doesn't
if [ ! -f "$CHR_SIZES" ]; then
    echo "Chromosome sizes file not found. Downloading..."
    bash $MODULE_PATH/download_chr_sizes.sh $CHR_SIZES $ASSEMBLY
fi

# Check and convert chromosome sizes to bed file if not present
if [ ! -f "${CHR_SIZES}.bed" ]; then
    echo "Converting chromosome sizes to bed file"
    awk '{print $1 "\t0\t" $2}' $CHR_SIZES > "${CHR_SIZES}.bed"
fi

# Peak type - narrow or broad use different configurations
if [ $USE_SINGULARITY == "True" ]; then
    echo "Using singularity"
    IMAGE_URI="docker://fooliu/macs2:version-2.2.7.1"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif

    if [ $PEAK_TYPE == "narrow" ]; then
        # Call peaks using MACS2 for narrow peaks
        singularity exec --bind $BIND_PATH $CONTAINER \
        macs2 callpeak -t $IP_SAMPLE -c $INPUT_SAMPLE \
                    -q 0.01 \
                    -n macs2 --outdir $OUTPUT_PATH
        cp $OUTPUT_PATH/macs2_peaks.narrowPeak $OUTPUT_PEAK

    elif [ $PEAK_TYPE == "broad" ]; then
        # Call peaks using MACS2 for broad peaks
        singularity exec --bind $BIND_PATH $CONTAINER \
        macs2 callpeak -t $IP_SAMPLE -c $INPUT_SAMPLE \
                    --broad \
                    -q 0.01 \
                    -n macs2 --outdir $OUTPUT_PATH
        cp $OUTPUT_PATH/macs2_peaks.broadPeak $OUTPUT_PEAK

    else
        echo "Please specify peak type: narrow or broad"
        exit 1
    fi

    # Further processing with bedtools
    IMAGE_URI="docker://staphb/bedtools:2.30.0"
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    CONTAINER=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    singularity exec --bind $BIND_PATH $CONTAINER \
    bedtools intersect -a $OUTPUT_PEAK -b "${CHR_SIZES}.bed" > ${OUTPUT_PEAK}.filtered
    singularity exec --bind $BIND_PATH $CONTAINER \
    bedtools sort -i ${OUTPUT_PEAK}.filtered -g $CHR_SIZES > ${OUTPUT_PEAK}.sorted
    singularity exec --bind $BIND_PATH $CONTAINER \
    bedtools complement -i ${OUTPUT_PEAK}.sorted -g $CHR_SIZES > $OUTPUT_NON_PEAK
else
    module load macs2/2.1.1
    if [ $PEAK_TYPE == "narrow" ]; then
        macs2 callpeak -t $IP_SAMPLE -c $INPUT_SAMPLE \
                    -q 0.01 \
                    -n macs2 --outdir $OUTPUT_PATH
        cp $OUTPUT_PATH/macs2_peaks.narrowPeak $OUTPUT_PEAK

    elif [ $PEAK_TYPE == "broad" ]; then
        macs2 callpeak -t $IP_SAMPLE -c $INPUT_SAMPLE \
                    --broad \
                    -q 0.01 \
                    -n macs2 --outdir $OUTPUT_PATH
        cp $OUTPUT_PATH/macs2_peaks.broadPeak $OUTPUT_PEAK

    else
        echo "Please specify peak type: narrow or broad"
        exit 1
    fi

    module load bedtools/2.30.0
    bedtools intersect -a $OUTPUT_PEAK -b "${CHR_SIZES}.bed" > ${OUTPUT_PEAK}.filtered
    bedtools sort -i ${OUTPUT_PEAK}.filtered -g $CHR_SIZES > ${OUTPUT_PEAK}.sorted
    bedtools complement -i ${OUTPUT_PEAK}.sorted -g $CHR_SIZES > $OUTPUT_NON_PEAK
fi
