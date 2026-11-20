#!/bin/bash

cd /COMMA/001_scsnRNA-seq/001_integration_benchmarking/001_integration_secondaries 

# Handy script to submit all secondary integration scripts together

for script in *.sh
do
echo "Submitting $script to cluster"
sbatch $script -r "/COMMA/results/scrna/integrated/pre_integration/preprocessed.h5ad" -h "/COMMA/results/scrna/integrated/pre_integration/preprocessed_hvg.h5ad"
done
