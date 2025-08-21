#!/bin/bash

peak_file=$3
non_peak_file=$4

parent_path=$(dirname $peak_file)
mkdir -p $parent_path

cd $parent_path
# Download files
# Check if the file exists
if [ ! -f cPeaks_hg38_raw.bed ]; then
    echo "Downloading cPeaks_hg38_raw.bed"
    #wget -O cPeaks_hg38_raw.bed https://github.com/MengQiuchen/cPeaks/raw/main/cPeaks_hg38.bed
    wget -O cPeaks_hg38_raw.bed https://cloud.tsinghua.edu.cn/f/6eb530748b324f53bc1f/?dl=1
fi

if [ ! -f hg38.chrom.sizes ]; then
    echo "Downloading hg38.chrom.sizes"
    wget -O hg38.chrom.sizes https://raw.githubusercontent.com/igvteam/igv/master/genomes/sizes/hg38.chrom.sizes
fi

function remove_excessive_chromosomes() {
    # Input file
    input_file=hg38.chrom.sizes

    # Output file
    output_file=hg38.autosomes.chrX.chrom.sizes

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

chrom_sizes=hg38.autosomes.chrX.chrom.sizes

# Complement peak file
if [ $1 == "True" ]; then
    # Sort the bed file
    singularity exec --bind $2  docker://staphb/bedtools:2.30.0 \
    grep -v 'chrY' cPeaks_hg38_raw.bed > cPeaks_hg38_no_chrY.bed 
    singularity exec --bind $2  docker://staphb/bedtools:2.30.0 \
    bedtools sort -i cPeaks_hg38_no_chrY.bed -g $chrom_sizes > $peak_file
    singularity exec --bind $2  docker://staphb/bedtools:2.30.0 \
    bedtools complement -i $peak_file -g $chrom_sizes > $non_peak_file
else
    module load bedtools/2.30.0
    grep -v 'chrY' cPeaks_hg38_raw.bed > cPeaks_hg38_no_chrY.bed 
    bedtools sort -i cPeaks_hg38_no_chrY.bed -g $chrom_sizes > $peak_file
    bedtools complement -i $peak_file -g $chrom_sizes > $non_peak_file
fi
