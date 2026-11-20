#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "scanpy environment loaded"

cd /COMMA/results/scrna/downstream_analysis/subclustering

python /COMMA/002_VSMC_phenotypes/004_PDGF_signature/005_cluster_synthetic_score.sh