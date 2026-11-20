#!/bin/bash

i=$((SLURM_ARRAY_TASK_ID-1))

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

InputFolders=$(ls -d /COMMA/data/meta_analysis/processed/*/*Solo.out/GeneFull_Ex50pAS/ | grep -v "spatial") # Exclude spatial samples from RNA meta-analysis

# Set up array to run in parallel
fnames=()
for fname in $InputFolders; do
    fnames+=("${fname%Solo.out/GeneFull_Ex50pAS/}")
done

cd ${fnames[$i]}Solo.out/

# Run scrublet for each sample & save to .h5ad output
echo "Submitting Python Job"
python /COMMA/000_preprocessing/003_RNA_QC/Python_scripts/002_run_scrublet.py -f ${fnames[$i]}Solo.out/GeneFull_Ex50pAS/cellbender/cellbender_output_filtered.h5
