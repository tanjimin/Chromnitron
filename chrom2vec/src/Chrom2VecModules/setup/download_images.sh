#!/bin/bash

# Function to check if a Singularity image exists locally
check_local_singularity_image() {
    IMAGE_PATH="$1"
    if [ -f "$IMAGE_PATH" ]; then
        echo "True"
    else
        echo "False"
    fi
}

# Function to pull Singularity image if not exists locally
pull_singularity_image_if_not_exists() {
    IMAGE_URI="$1"
    # Remove 'docker://' prefix if exists and replace all '/' with '-'
    LOCAL_IMAGE_PATH="${IMAGE_URI/docker:\/\//}"
    LOCAL_IMAGE_PATH=~/.singularity/cache/"${LOCAL_IMAGE_PATH//\//_}".sif
    mkdir -p ~/.singularity/cache/
    if [ $(check_local_singularity_image "$LOCAL_IMAGE_PATH") == "False" ]; then
        echo "Pulling image: $IMAGE_URI"
        singularity pull "$LOCAL_IMAGE_PATH" "$IMAGE_URI"
    else
        echo "Using local image: $LOCAL_IMAGE_PATH"
        
    fi
    echo "$LOCAL_IMAGE_PATH"
}

# Pull all the images
pull_singularity_image_if_not_exists "docker://staphb/samtools:1.17" 
pull_singularity_image_if_not_exists "docker://biocontainers/hisat2:v2.1.0-2-deb_cv1"
pull_singularity_image_if_not_exists "docker://stjudecloud/deeptools:branch-encoding-1.0.0"
pull_singularity_image_if_not_exists "docker://juettemann/genrich:0.6"
pull_singularity_image_if_not_exists "docker://tanjimin/seq2vec:0.0.1"
pull_singularity_image_if_not_exists "docker://tanjimin/pybigwig:0.0.1"
pull_singularity_image_if_not_exists "docker://fooliu/macs2:version-2.2.7.1"
pull_singularity_image_if_not_exists "docker://staphb/bedtools:2.30.0"
pull_singularity_image_if_not_exists "docker://staphb/fastqc:0.12.1"
pull_singularity_image_if_not_exists "docker://staphb/fastp:0.23.4"
