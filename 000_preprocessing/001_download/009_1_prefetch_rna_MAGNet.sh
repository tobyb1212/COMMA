#!/bin/bash

# Downloading Files from dbGaP phs001539.v4.p1

cd /COMMA/data/meta_analysis/fastq/rna_magnet

export PATH="$PATH:software/sratoolkit.3.1.1-centos_linux64/bin"
while read -r file; do
  prefetch $file --ngc magnet.ngc --max-size 100g
done < accession_list.txt
