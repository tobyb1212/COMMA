#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scarches
echo "scanpy environment loaded"

cd /COMMA/results/scrna/downstream_analysis/subclustering

python /COMMA/001_scsnRNA-seq/003_vascular_integration/001_vascular_reintegrate_secondaries/Python_scripts/001_vascular_integrate_scPoli.py