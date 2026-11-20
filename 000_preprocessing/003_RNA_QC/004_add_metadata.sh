#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

echo "Submitting Python Job"
python /COMMA/000_preprocessing/003_RNA_QC/Python_scripts/004_add_metadata.py