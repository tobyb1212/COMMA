#!/bin/bash

module load Miniforge3
eval "$(conda shell.bash hook)"
conda activate snakemake 
module load R

cd software/pyflow_seurat_parameter_adapted

sbatch software/pyflow_seurat_parameter_adapted/pyflow-scBoot.sh

echo "Plotting R"
Rscript /COMMA/002_VSMC_phenotypes/002_cluster_stability/R_scripts/001_cluster_stability.R