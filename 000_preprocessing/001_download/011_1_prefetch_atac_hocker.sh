#!/bin/bash
# Download SRA Files from phs002204.v1.p1

cd /COMMA/data/meta_analysis/fastq/atac_hocker

export PATH="$PATH:software/sratoolkit.3.1.1-centos_linux64/bin"
while read -r file; do
  prefetch $file --ngc hocker.ngc --max-size 100g 
done < accession_list.txt