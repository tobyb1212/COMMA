#!/bin/bash

export PATH="$PATH:software/HiCUP-0.9.2"
export PATH="$PATH:software/bowtie2" # Artificially Setting the Path as Bowtie Bin
module load R-bundle-Bioconductor/3.18-foss-2023a-R-4.3.2
module load SAMtools/1.18-GCC-12.3.0

hicup --config /COMMA/results/CCC/CA-VSMCs/scripts/HiCUP_CASMC.txt