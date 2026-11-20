#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scvi-env
echo "scanpy environment loaded"

cd /COMMA/results/scrna/downstream_analysis/subclustering/vsmc

python /COMMA/002_VSMC_phenotypes/001_VSMC_integration/Python_scripts/001_vsmc_integrate_scVI.py