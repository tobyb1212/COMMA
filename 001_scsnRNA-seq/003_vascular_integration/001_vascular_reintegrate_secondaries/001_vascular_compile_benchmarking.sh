#!/bin/bash

cd /COMMA/results/scrna/downstream_analysis/subclustering/vascular/
module load Miniforge3

eval "$(conda shell.bash hook)"
conda activate scanpy
echo "Scanpy environment loaded"

InputFiles=$(ls *_metrics.csv | tr '\n' ',')
InputFiles=${InputFiles%,}

echo "Submitting Python Job"
python /COMMA/001_scsnRNA-seq/003_vascular_integration/001_vascular_reintegrate_secondaries/Python_scripts/001_vascular_compile_benchmarking.py -f $InputFiles