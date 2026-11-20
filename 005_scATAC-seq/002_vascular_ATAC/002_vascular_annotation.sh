#!/bin/bash

module load R 
module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

Rscript /COMMA/005_scATAC-seq/002_vascular_ATAC/R_scripts/002_vascular_annotation.R

python /COMMA/005_scATAC-seq/002_vascular_ATAC/Python_scripts/002_vascular_annotation.py