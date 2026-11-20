#!/bin/bash

module load R 

InputFiles=$(ls -d /COMMA/data/meta_analysis/processed/*/*/outs/*ragments.tsv.gz)

Rscript /COMMA/000_preprocessing/004_ATAC_QC/R_scripts/001_create_ArchR.R -f "$InputFiles"