#!/bin/bash

cd /COMMA/results/scrna/integrated
module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

InputFiles=$(ls */*_metrics.csv | tr '\n' ',')
InputFiles=${InputFiles%,}

echo "Submitting Python Job"
python /COMMA/001_scsnRNA-seq/001_integration_benchmarking/Python_scripts/005_compile_integration_benchmarking.py -f $InputFiles