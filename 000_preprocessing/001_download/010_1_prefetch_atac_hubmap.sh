#!/bin/bash
# Downloading Files from dbGaP phs002267.v2.p1 

cd /COMMA/data/meta_analysis/fastq/atac_hubmap

export PATH="$PATH:software/sratoolkit.3.1.1-centos_linux64/bin"
while read -r file; do
  prefetch $file --ngc hubmap.ngc --max-size 100g
done < "accession_list_new.txt"