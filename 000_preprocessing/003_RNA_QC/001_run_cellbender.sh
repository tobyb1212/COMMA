#!/bin/bash

i=$((SLURM_ARRAY_TASK_ID-1))

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

InputFolders=$(ls -d /COMMA/data/meta_analysis/processed/*/*Solo.out/GeneFull_Ex50pAS/ | grep -v "spatial") # Exclude spatial samples from RNA meta-analysis

# Set up array job so cellbender runs in parallel
fnames=()
for fname in $InputFolders; do
    fnames+=("${fname%Solo.out/GeneFull_Ex50pAS/}")
done

cd ${fnames[$i]}Solo.out/

# Create .h5ad Object from gzipped STARsolo output
echo "Submitting Python Job"
python /COMMA/000_preprocessing/003_RNA_QC/Python_scripts/001_run_cellbender.py -f ${fnames[$i]}Solo.out/GeneFull_Ex50pAS/raw_gzip/

conda deactivate 
conda activate cellbender

# Run CellBender on Each Individual Sample
cellbender remove-background \
     --input ${fnames[$i]}Solo.out/GeneFull_Ex50pAS/cellbender/raw_feature_barcode_matrix.h5ad \
     --output ${fnames[$i]}Solo.out/GeneFull_Ex50pAS/cellbender/cellbender_output.h5