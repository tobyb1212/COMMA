#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "scanpy environment loaded"

python /COMMA/004_CVD_DGE/Python_scripts/001_dge_analysis.py