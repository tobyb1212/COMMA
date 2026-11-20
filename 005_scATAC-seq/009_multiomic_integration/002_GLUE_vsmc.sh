#!/bin/bash

module load R
module load Miniforge3
eval "$(conda shell.bash hook)"

conda activate GLUE

# Multiomic Integration Analysis
python /COMMA/005_scATAC-seq/009_multiomic_integration/Python_scripts/002_GLUE_vsmc.py

# TF Inference
python /COMMA/005_scATAC-seq/009_multiomic_integration/Python_scripts/002_TF_inference.py