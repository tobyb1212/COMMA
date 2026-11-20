#!/bin/bash

module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scvi-env
echo "Scanpy environment loaded"

cd "/COMMA/results/scrna/downstream_analysis/epicardial_coronary/integrated" || exit

echo "Submitting Python Job"
python //COMMA/003_Epicardial_Comparison/002_integration/Python_scripts/001_integration.py