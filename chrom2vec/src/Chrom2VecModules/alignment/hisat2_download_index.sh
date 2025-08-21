#!/bin/bash

# Assigning positional parameters to descriptive variable names
USE_SINGULARITY=$1
BIND_PATH=$2
INDEX_ROOT=$3
ASSEMBLY=$4

mkdir -p $INDEX_ROOT 
cd $INDEX_ROOT
wget "https://genome-idx.s3.amazonaws.com/hisat/${ASSEMBLY}_genome.tar.gz"
tar -xvzf "${ASSEMBLY}_genome.tar.gz"
