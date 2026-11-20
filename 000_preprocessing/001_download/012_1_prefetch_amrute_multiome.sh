#!/bin/bash
# Downloading SRA Files from SRP516110

cd /COMMA/data/meta_analysis/fastq/multiome_amrute/

export PATH="$PATH:software/sratoolkit.3.1.1-centos_linux64/bin"
while read -r file; do
  prefetch $file --max-size 100g 
done < /COMMA/data/meta_analysis/fastq/multiome_amrute/failed_starsolo.txt