#!/bin/bash

i=$((SLURM_ARRAY_TASK_ID-1))

cd /COMMA/data/meta_analysis/fastq/multiome_amrute/

export PATH="$PATH:software/sratoolkit.3.1.1-centos_linux64/bin"

# Extracting SRA Files from SRA SRP516110
# Including Techincal Reads for ATAC-seq

fnames=()
for fname in */
do
fnames[${#fnames[@]}]=${fname%/}
done

fasterq-dump ${fnames[$i]} --outdir /COMMA/data/meta_analysis/fastq/multiome_amrute/ --split-files --include-technical
    if [[ $? -ne 0 ]]; then
        echo "Error: fasterq-dump failed for ${fnames[$i]}"
    fi

cd ${fnames[$i]}
rm ${fnames[$i]}.sra
cd /COMMA/data/meta_analysis/fastq/multiome_amrute/
gzip ${fnames[$i]}*.fastq -4

cd /COMMA/data/meta_analysis/fastq/multiome_amrute/

# Renaming Files According to CellRanger Naming Conventions
# Change for RNA vs ATAC

for gz_file in *.fastq.gz; do
        if [[ "$gz_file" == *_1.fastq.gz ]]; then
            mv "$gz_file" "${gz_file/_1/_S1_L001_I1_001}"
        elif [[ "$gz_file" == *_2.fastq.gz ]]; then
            mv "$gz_file" "${gz_file/_2/_S1_L001_R1_001}"
        elif [[ "$gz_file" == *_3.fastq.gz ]]; then
            mv "$gz_file" "${gz_file/_3/_S1_L001_R2_001}"
        fi
    done

