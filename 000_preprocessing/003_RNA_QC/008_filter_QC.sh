#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

cd "/COMMA/results/scrna/count_matrix" || exit

echo "Submitting Python Job"
# Submitting anndata object to be filtered for QC
python /COMMA/000_preprocessing/003_RNA_QC/Python_scripts/008_filter_QC.py

