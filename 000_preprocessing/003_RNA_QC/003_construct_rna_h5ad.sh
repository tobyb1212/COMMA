#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

InputFiles=$(ls /COMMA/data/meta_analysis/processed/*/*Solo.out/GeneFull_Ex50pAS/cellbender/scrublet.h5ad | grep -v "spatial" | tr '\n' ',') # Excluding Spatial Samples from snRNA-seq Meta-Analysis
OutDir=/COMMA/results/scrna/count_matrix
Summary_csvs=$(ls /COMMA/data/meta_analysis/processed/*/*Solo.out/GeneFull_Ex50pAS/Summary.csv | grep -v "spatial" | tr '\n' ',') # Submitting Alignment Summary Statistics to Command Line

# Drop any trailing commas
InputFiles=${InputFiles%,}
Summary_csvs=${Summary_csvs%,}

cd $OutDir || exit

echo "Submitting Python Job"
python /COMMA/000_preprocessing/003_RNA_QC/Python_scripts/003_construct_rna_h5ad.py -f $InputFiles -o $OutDir -s $Summary_csvs