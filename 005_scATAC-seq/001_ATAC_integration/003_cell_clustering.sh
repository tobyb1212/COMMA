#!/bin/bash

module load R 
module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

Rscript /COMMA/005_scATAC-seq/001_ATAC_integration/R_scripts/003_cluster_annotations.R

python /COMMA/005_scATAC-seq/001_ATAC_integration/Python_scripts/003_cell_clustering.py