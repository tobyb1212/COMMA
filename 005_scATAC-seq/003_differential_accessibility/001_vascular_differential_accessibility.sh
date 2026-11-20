#!/bin/bash

module load R
module load Miniforge3
eval "$(conda shell.bash hook)"

# Submitting R Script to Convert 
Rscript /COMMA/005_scATAC-seq/003_differential_accessibility/R_scripts/001_vascular_differential_accessibility.R
