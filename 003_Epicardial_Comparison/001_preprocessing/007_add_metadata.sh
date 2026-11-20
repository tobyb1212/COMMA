#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

echo "Submitting Python Job"
python /COMMA/003_Epicardial_Comparison/001_preprocessing/Python_scripts/007_add_metadata.py