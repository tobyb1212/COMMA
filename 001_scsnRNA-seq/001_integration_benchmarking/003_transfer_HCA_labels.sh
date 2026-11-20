#!/bin/bash

i=$((SLURM_ARRAY_TASK_ID-1))

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

cd /COMMA/results/scrna/integrated/

# Generate an array to label cells in parallel
fnames=()
for fname in */*_doublet_removed.h5ad; do
    fnames+=("${fname%_doublet_removed.h5ad}")
done

# For each integration methodology cluster and remove doublets
echo "Submitting Python Job"
python /COMMA/001_scsnRNA-seq/001_integration_benchmarking/Python_scripts/003_transfer_HCA_labels.py -f ${fnames[$i]}_doublet_removed.h5ad