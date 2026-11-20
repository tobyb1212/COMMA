#!/bin/bash

i=$((SLURM_ARRAY_TASK_ID-1))

module load Miniforge3
module load R

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "scanpy environment loaded"

cd /COMMA/results/scrna/downstream_analysis/subclustering/vascular

# Generate an array to label cells in parallel
fnames=()
for fname in *_vascular_integrated.h5ad; do
    fnames+=("${fname%_vascular_integrated.h5ad}")
done

rawfile="/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular.h5ad"
integrated_file="/COMMA/results/scrna/downstream_analysis/subclustering/vascular/${fnames[$i]}_vascular_integrated.h5ad"

python /COMMA/001_scsnRNA-seq/003_vascular_integration/001_vascular_reintegrate_secondaries/Python_scripts/001_vascular_integration_evaluation.py -i $integrated_file -f $rawfile
