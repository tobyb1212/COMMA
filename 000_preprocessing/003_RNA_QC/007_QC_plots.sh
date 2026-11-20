#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

# Finding individual study count matrices
InputFiles=$(ls /COMMA/results/scrna/count_matrix/study_matrices/*.h5ad | tr '\n' ',') 
InputFiles=${InputFiles%,} # Drop any trailing commas

cd /COMMA/results/scrna/QC || exit

# Submitting python script to generate QC plots
python /COMMA/000_preprocessing/003_RNA_QC/Python_scripts/007_QC_plots.py -f $InputFiles