#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

cd "/COMMA/results/scrna/count_matrix" || exit

echo "Submitting Python Job"
# Submitting anndata object for processing before integration
python /COMMA/000_preprocessing/003_RNA_QC/Python_scripts/009_preprocess.py

