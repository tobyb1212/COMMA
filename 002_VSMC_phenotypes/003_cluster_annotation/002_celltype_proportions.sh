#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

echo "Submitting Python Job"
python /COMMA/002_VSMC_phenotypes/003_cluster_annotation/Python_scripts/002_celltype_proportions.py