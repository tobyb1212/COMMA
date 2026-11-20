#!/bin/bash

i=$((SLURM_ARRAY_TASK_ID-1))

module load Miniforge3
module load R

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "scanpy environment loaded"

cd /COMMA/results/scrna/integrated/

# Generate an array to label cells in parallel
fnames=()
for fname in */*_celltypist_predicted.h5ad; do
    fnames+=("${fname%_celltypist_predicted.h5ad}")
done

rawfile="/COMMA/results/scrna/integrated/pre_integration/preprocessed.h5ad"

python /COMMA/001_scsnRNA-seq/001_integration_benchmarking/Python_scripts/004_benchmark_integration.py -i ${fnames[$i]}_celltypist_predicted.h5ad -f $rawfile
