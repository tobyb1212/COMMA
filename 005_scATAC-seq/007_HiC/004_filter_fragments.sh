#!/bin/bash

module load R
module load SAMtools
module load HTSlib

Rscript /COMMA/005_scATAC-seq/007_HiC/R_scripts/004_filter_fragments_files.R