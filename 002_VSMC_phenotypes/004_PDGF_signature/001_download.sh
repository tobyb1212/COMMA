#!/bin/bash

cd /COMMA/data/RNA/pdgfbb/conley

export PATH="$PATH:software/sratoolkit.3.1.1-centos_linux64/bin"
while read -r file; do
  prefetch $file --max-size 100g 
done < accession_list.txt

for file in */
do 

fasterq-dump $file --outdir /COMMA/data/RNA/pdgfbb/conley --split-files --include-technical
    if [[ $? -ne 0 ]]; then
        echo "Error: fasterq-dump failed for $file"
        continue
    fi
cd $file
name=$(basename $file)
rm "${name}.sra"
cd /COMMA/data/RNA/pdgfbb/conley
gzip ${name}*.fastq -4

done

cd /COMMA/data/RNA/pdgfbb/conley


