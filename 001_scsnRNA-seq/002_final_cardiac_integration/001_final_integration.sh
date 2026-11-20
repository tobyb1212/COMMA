#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

cd /COMMA/results/scrna/integrated/

echo "Submitting Python Job"
python /COMMA/001_scsnRNA-seq/002_final_cardiac_integration/Python_scripts/001_final_integration.py