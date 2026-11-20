#!/bin/bash

module load Miniforge3
eval "$(conda shell.bash hook)"
module load R

Rscript /COMMA/005_scATAC-seq/001_ATAC_integration/R_scripts/001_transfer_multiome_labels.R