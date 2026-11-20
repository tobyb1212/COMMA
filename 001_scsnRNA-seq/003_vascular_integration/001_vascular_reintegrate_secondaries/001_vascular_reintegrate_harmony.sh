#!/bin/bash

module load Miniforge3
module load R

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "scanpy environment loaded"

cd /COMMA/results/scrna/downstream_analysis/subclustering/vascular

Rscript /COMMA/001_scsnRNA-seq/003_vascular_integration/001_vascular_reintegrate_secondaries/R_scripts/001_vascular_integrate_harmony.R

python /COMMA/001_scsnRNA-seq/003_vascular_integration/001_vascular_reintegrate_secondaries/Python_scripts/001_vascular_integrate_harmony.py