#!/bin/bash

while getopts "r:h:" option; do 
   case "$option" in 
        r) RawFile=${OPTARG} ;;
        h) hvgFile=${OPTARG} ;;
        *) echo "Usage: $0 [-r RawFile] [-h HVG File]" >&2
        exit 1 ;;
    esac
done

if [ -f "$RawFile" ]
then
        echo "Raw file $RawFile found"
else
        echo "Raw file not found"
        exit 0
fi

if [ -f "$hvgFile" ]
then
        echo "Raw file $hvgFile found"
else
        echo "Raw file not found"
        exit 0
fi

module load R
module load Miniforge3
eval "$(conda shell.bash hook)"

# Submitting Non-Integrated Count Matrix to R CCA Integration Script & Python Evaluation Script
Rscript /COMMA/001_scsnRNA-seq/001_integration_benchmarking/001_integration_secondaries/R_scripts/001_06_integration_harmony.R -v $hvgFile

conda activate scvi-env

python /COMMA/001_scsnRNA-seq/001_integration_benchmarking/001_integration_secondaries/Python_scripts/001_06_integration_harmony.py -f $RawFile