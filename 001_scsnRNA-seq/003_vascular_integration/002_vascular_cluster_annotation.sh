#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "scanpy environment loaded"

cd /COMMA/results/scrna/downstream_analysis/subclustering/vascular

python /COMMA/001_scsnRNA-seq/003_vascular_integration/Python_scripts/002_vascular_cluster_annotation.py