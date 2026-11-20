#!/bin/bash
# Downloading SRA Files from SRA SRP422510

cd /COMMA/data/meta_analysis/fastq/atac_selewa/

export PATH="$PATH:software/sratoolkit.3.1.1-centos_linux64/bin"
while read -r file; do
  prefetch $file --max-size 100g 
done < /COMMA/data/meta_analysis/fastq/atac_selewa/accession_list.txt