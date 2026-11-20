#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

cd /COMMA/results/scrna/downstream_analysis/epicardial_coronary/count_matrix/

file="/COMMA/results/scrna/downstream_analysis/epicardial_coronary/count_matrix/rna_count_matrix_metadata.h5ad"

python /COMMA/003_Epicardial_Comparison/001_preprocessing/Python_scripts/008_assign_genescores.py -f $file