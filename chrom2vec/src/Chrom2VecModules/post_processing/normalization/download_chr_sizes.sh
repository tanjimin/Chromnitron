#!/bin/bash

chr_file=$1
assembly=$2
chr_sizes_raw=$assembly.chrom.sizes

parent_path=$(dirname $chr_file)
mkdir -p $parent_path

cd $parent_path

#if [ ! -f hg38.chrom.sizes ]; then
if [ ! -f $chr_sizes_raw ]; then
    echo "Downloading $chr_sizes_raw"
    wget -O $chr_sizes_raw https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/$chr_sizes_raw
fi

function remove_excessive_chromosomes() {
    # Input file
    input_file=$chr_sizes_raw

    # Output file
    output_file=$chr_file

    # Empty the output file if it already exists
    > $output_file

    # Iterate over each line in the file
    while IFS= read -r line
    do
        # Extract the chromosome name (assuming it's the first field in the line)
        chrom=$(echo "$line" | cut -f1)

        # If the chromosome name matches autosomes or chrX, append it to the output file
        if [[ "$chrom" =~ ^chr([1-9]|1[0-9]|2[0-2]|X)$ ]]; then
            echo "$line" >> "$output_file"
        fi
    done < "$input_file"
}

remove_excessive_chromosomes
