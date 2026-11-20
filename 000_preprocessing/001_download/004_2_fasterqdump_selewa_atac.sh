#!/bin/bash

cd /COMMA/data/meta_analysis/fastq/atac_selewa/

export PATH="$PATH:software/sratoolkit.3.1.1-centos_linux64/bin"

# Extracting all Files (Including Technical Index Reads for ATAC-seq) for PRJNA934593
# Files are gzipped for storage

for file in */
do 

fasterq-dump $file --outdir /COMMA/data/meta_analysis/fastq/atac_selewa/ --split-files --include-technical
    if [[ $? -ne 0 ]]; then
        echo "Error: fasterq-dump failed for $file"
        continue
    fi
cd $file
name=$(basename $file)
rm "${name}.sra"
cd /COMMA/data/meta_analysis/fastq/atac_selewa/
gzip ${name}*.fastq -4

done

cd /COMMA/data/meta_analysis/fastq/atac_selewa/

# Rename Files to Reflect CellRanger Naming Conventions

for gz_file in *.fastq.gz; do
        if [[ "$gz_file" == *_1.fastq.gz ]]; then
            mv "$gz_file" "${gz_file/_1/_S1_L001_I1_001}"
        elif [[ "$gz_file" == *_2.fastq.gz ]]; then
            mv "$gz_file" "${gz_file/_2/_S1_L001_R1_001}"
        elif [[ "$gz_file" == *_3.fastq.gz ]]; then
            mv "$gz_file" "${gz_file/_3/_S1_L001_I2_001}"
        elif [[ "$gz_file" == *_4.fastq.gz ]]; then
            mv "$gz_file" "${gz_file/_4/_S1_L001_R2_001}"
        fi
    done

