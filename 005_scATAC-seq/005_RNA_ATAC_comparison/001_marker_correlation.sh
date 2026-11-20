#!/bin/bash

module load R
module load Miniforge3
eval "$(conda shell.bash hook)"

# Submitting R Script to Convert 
Rscript /COMMA/005_scATAC-seq/005_RNA_ATAC_comparison/R_scripts/001_marker_correlation.R

conda activate scanpy

python /COMMA/005_scATAC-seq/005_RNA_ATAC_comparison/Python_scripts/001_marker_correlation.py