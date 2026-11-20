#!/bin/bash

module load R
module load Miniforge3
eval "$(conda shell.bash hook)"
conda activate scanpy

# Submitting R Script to Convert 
Rscript /COMMA/005_scATAC-seq/002_vascular_ATAC/R_scripts/001_vascular_reintegrate.R

conda deactivate
conda activate scvi-env

python /COMMA/005_scATAC-seq/002_vascular_ATAC/Python_scripts/001_vascular_PeakVI.py