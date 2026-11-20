#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

cd "/COMMA/results/scrna/count_matrix" || exit

# Splitting adata object by study to make processing easier
echo "Submitting Python Job"
python /COMMA/000_preprocessing/003_RNA_QC/Python_scripts/005_split_adata.py