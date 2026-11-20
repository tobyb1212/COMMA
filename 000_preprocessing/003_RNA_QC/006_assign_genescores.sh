#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

cd /COMMA/results/scrna/count_matrix/study_matrices/

# Add Stress & Cell Cycle Gene Scores for Each Individual Split .h5ad Object
for file in *.h5ad; do
python /COMMA/000_preprocessing/003_RNA_QC/Python_scripts/006_assign_gene_scores.py -f $file
done