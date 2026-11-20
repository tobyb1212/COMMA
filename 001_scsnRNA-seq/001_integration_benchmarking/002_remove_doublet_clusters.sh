#!/bin/bash

i=$((SLURM_ARRAY_TASK_ID-1))

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

cd /COMMA/results/scrna/integrated/

# Generate an array to remove doublets in parallel
fnames=()
for fname in */*integrated.h5ad; do
    fnames+=("${fname%.h5ad}")
done

# For each integration methodology cluster and remove doublets
echo "Submitting Python Job"
python /COMMA/001_scsnRNA-seq/001_integration_benchmarking/Python_scripts/002_remove_doublet_clusters.py -f ${fnames[$i]}.h5ad
