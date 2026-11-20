#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "scanpy environment loaded"

cd /COMMA/results/scrna/downstream_analysis/subclustering

python /COMMA/004_CVD_DGE/Python_scripts/002_pseudobulk_dge.py

module purge
module load R

Rscript /COMMA/004_CVD_DGE/R_scripts/002_pseudobulk_dge.R
Rscript /COMMA/004_CVD_DGE/R_scripts/002_consensus_dge.R