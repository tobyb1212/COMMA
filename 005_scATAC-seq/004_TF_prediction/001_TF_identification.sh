#!/bin/bash

module load R
module load GSL
module load Miniforge3
eval "$(conda shell.bash hook)"

Rscript /COMMA/005_scATAC-seq/004_TF_prediction/R_scripts/001_TF_identification.R

conda activate scanpy
python /COMMA/005_scATAC-seq/004_TF_prediction/Python_scripts/001_sox_dotplot.py
