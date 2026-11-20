#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

cd "/COMMA/results/scrna/downstream_analysis/epicardial_coronary/DGE" || exit

echo "Submitting Python Job"    
python /COMMA/003_Epicardial_Comparison/004_DGE/python_scripts/001_vsmc_dge.py

module load R

Rscript /COMMA/003_Epicardial_Comparison/004_DGE/R_scripts/001_vsmc_dge.R
Rscript /COMMA/003_Epicardial_Comparison/004_DGE/R_scripts/001_consensus_dge.R
Rscript /COMMA/003_Epicardial_Comparison/004_DGE/R_scripts/001_vsmc_pathway.R
