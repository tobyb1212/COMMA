#!/bin/bash

cd /COMMA/data/meta_analysis/fastq/atac_hubmap

export PATH="$PATH:software/sratoolkit.3.1.1-centos_linux64/bin"

# Extracting all Reads from SRA Files of HUBMAP
# Using ngc Code for dbGaP Download

for file in SRR*/
do 

fasterq-dump $file --outdir /COMMA/data/meta_analysis/fastq/atac_hubmap --ngc hubmap.ngc --split-files --include-technical
    if [[ $? -ne 0 ]]; then
        echo "Error: fasterq-dump failed for $file"
        continue
    fi
cd $file
name=$(basename $file)
rm "${name}.sra"
cd /COMMA/data/meta_analysis/fastq/atac_hubmap
gzip ${name}*.fastq -4

done

cd /COMMA/data/meta_analysis/fastq/atac_hubmap

# Renaming Files to Fit Naming Conventions

for gz_file in *54*.fastq.gz; do
        if [[ "$gz_file" == *_1.fastq.gz ]]; then
            mv "$gz_file" "${gz_file/_1/_S1_L001_R1_001}"
        elif [[ "$gz_file" == *_2.fastq.gz ]]; then
            mv "$gz_file" "${gz_file/_2/_S1_L001_R2_001}"
        fi
done

