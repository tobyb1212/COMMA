#!/bin/bash

export PATH="$PATH:software/HiCUP-0.9.2"

hicup_digester --genome GRCh38 --re1 ^GATC,MboI \
    /COMMA/reference_genome/GENCODE/GRCh38.p14.GENCODE_release44/GRCh38.primary_assembly.genome.fa \
    --outdir /COMMA/reference_genome/bowtie/HiCUP_restriction_MboI