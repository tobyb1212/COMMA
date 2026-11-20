#!/bin/bash

cd /COMMA/data/meta_analysis/fastq/rna_magnet

export PATH="$PATH:software/sratoolkit.3.1.1-centos_linux64/bin"

# Extracting fastq Files from SRA & Gzipping for Storage
# Using ngc Code for dbGaP Download

for file in SRR*/
do 

fasterq-dump $file --outdir /COMMA/data/meta_analysis/fastq/rna_magnet --ngc magnet.ngc --split-files --include-technical
    if [[ $? -ne 0 ]]; then
        echo "Error: fasterq-dump failed for $file"
        echo "$file" >> "failed_fasterqdump.txt"
        exit 1
    fi

name=$(basename $file)
rm ${name}_1.fastq # Removing I1 Read - not necessary
gzip ${name}*.fastq -4

done

# Changing File Names to Fit CellRanger Naming Conventions

for gz_file in *.fastq.gz; do
        if [[ "$gz_file" == *_1.fastq.gz ]]; then
            mv "$gz_file" "${gz_file/_1/_S1_L001_R1_001}"
        elif [[ "$gz_file" == *_2.fastq.gz ]]; then
            mv "$gz_file" "${gz_file/_2/_S1_L001_R2_001}"
        fi
    done

