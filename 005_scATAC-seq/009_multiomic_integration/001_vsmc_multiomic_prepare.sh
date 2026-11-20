#!/bin/bash

cd /COMMA/results/scatac/multiomic_integration/vsmc

module load Miniforge3
module load R
eval "$(conda shell.bash hook)"
conda activate scanpy

Rscript /COMMA/005_scATAC-seq/009_multiomic_integration/R_scripts/001_vsmc_multiomic_prepare.R