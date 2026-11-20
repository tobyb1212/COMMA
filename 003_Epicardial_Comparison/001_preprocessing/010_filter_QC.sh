#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

cd /COMMA/results/scrna/downstream_analysis/epicardial_coronary/count_matrix/

echo "Submitting Python Job"
python /COMMA/003_Epicardial_Comparison/001_preprocessing/Python_scripts/010_filter_QC.py
