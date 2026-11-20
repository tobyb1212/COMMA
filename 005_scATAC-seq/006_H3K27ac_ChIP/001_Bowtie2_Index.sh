#!/bin/bash

cd /COMMA/results/chIP/reference_genome || exit
module load Bowtie2/
bowtie2-build /COMMA/reference_genome/GENCODE/GRCh38.p14.GENCODE_release46/GRCh38.primary_assembly.genome.fa hg38 -p 2 