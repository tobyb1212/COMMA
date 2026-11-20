#!/bin/bash

module load R
module load Miniforge3
eval "$(conda shell.bash hook)"
conda activate scanpy

# Submitting R Script to Call Peaks & Convert to .H5AD
Rscript /COMMA/005_scATAC-seq/001_ATAC_integration/R_scripts/002_peakVI_final.R

conda deactivate
conda activate scvi-env

python COMMA/005_scATAC-seq/001_ATAC_integration/Python_scripts/002_peakVI.py