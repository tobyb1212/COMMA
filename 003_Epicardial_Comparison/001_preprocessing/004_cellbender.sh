#!/bin/bash

i=$((SLURM_ARRAY_TASK_ID-1))

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

InputFolders=$(ls -d /COMMA/results/scrna/downstream_analysis/epicardial_coronary/processed/*Solo.out/GeneFull_Ex50pAS/)

# Set up array job so cellbender runs in parallel
fnames=()
for fname in $InputFolders; do
    fnames+=("${fname%Solo.out/GeneFull_Ex50pAS/}")
done

cd ${fnames[$i]}Solo.out/

d=GeneFull_Ex50pAS/raw
cp -r "$d" "${d%/raw}/raw_gzip" && find "${d%/raw}/raw_gzip" -type f -exec gzip {} \;

# Create .h5ad Object from gzipped STARsolo output
echo "Submitting Python Job"
python /COMMA/003_Epicardial_Comparison/001_preprocessing/Python_scripts/004_run_cellbender.py -f ${fnames[$i]}Solo.out/GeneFull_Ex50pAS/raw_gzip/

conda deactivate 
conda activate cellbender

# Run CellBender on Each Individual Sample
cellbender remove-background \
     --input ${fnames[$i]}Solo.out/GeneFull_Ex50pAS/cellbender/raw_feature_barcode_matrix.h5ad \
     --output ${fnames[$i]}Solo.out/GeneFull_Ex50pAS/cellbender/cellbender_output.h5