#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

InputFiles=$(ls /COMMA/results/scrna/downstream_analysis/epicardial_coronary/processed/*Solo.out/GeneFull_Ex50pAS/cellbender/scrublet.h5ad | tr '\n' ',') # Excluding Spatial Samples from snRNA-seq Meta-Analysis
OutDir=/COMMA/results/scrna/downstream_analysis/epicardial_coronary/count_matrix
Summary_csvs=$(ls /COMMA/results/scrna/downstream_analysis/epicardial_coronary/processed/*Solo.out/GeneFull_Ex50pAS/Summary.csv | tr '\n' ',') # Submitting Alignment Summary Statistics to Command Line

# Drop any trailing commas
InputFiles=${InputFiles%,}
Summary_csvs=${Summary_csvs%,}

cd $OutDir || exit

echo "Submitting Python Job"
python /COMMA/003_Epicardial_Comparison/001_preprocessing/Python_scripts/006_construct_rna_h5ad.py -f $InputFiles -o $OutDir -s $Summary_csvs