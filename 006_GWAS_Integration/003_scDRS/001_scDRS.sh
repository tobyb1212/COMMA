#!/bin/bash

export PATH="software/magma_software:$PATH"

cd /COMMA/results/scrna/downstream_analysis/subclustering/vascular/scDRS

IN=/COMMA/results/GWAS/s-LDSC/GCST90132314_buildGRCh38.tsv

awk 'BEGIN{OFS="\t"; print "SNP","CHR","BP","P","N"}
NR>1 {
    print $23":"$24, $23, $24, $1, $21
}' "$IN" > GCST90132314_magma.tsv

# Collapse multiple genes of same symbol from GTF file (keep first chromosome)
awk 'BEGIN{OFS="\t"} $3=="gene" {
  match($0,/gene_name "[^"]+"/); g=substr($0,RSTART+11,RLENGTH-12); c=$1; sub(/^chr/,"",c)
  if(g!="" && c~/^([0-9]+|X|Y)$/ && !(g in seen)){seen[g]=1; print g,c,$4,$5,$7}
}' /COMMA/reference_genome/GENCODE/GRCh38.p14.GENCODE_release46/gencode.v46.basic.annotation_filtered.gtf \
> GRCh38_symbols.gene.loc

magma \
  --annotate \
  --snp-loc GCST90132314_magma.tsv \
  --gene-loc GRCh38_symbols.gene.loc \
  --out GCST90132314

magma \
    --bfile "/COMMA/reference_genome/ldsc/plink_map/merged/merged_ALL_chr" \
    --pval GCST90132314_magma.tsv ncol=N \
    --gene-annot GCST90132314.genes.annot \
    --out GCST90132314

awk 'BEGIN{OFS="\t"; print "GENE","CAD"}
NR>1 && $8!="NA" {
    print $1,$8
}' GCST90132314.genes.out > GCST90132314_gene.tsv

module load Miniforge3
eval "$(conda shell.bash hook)"
conda activate scdrs

scdrs munge-gs \
    --out-file /COMMA/results/scrna/downstream_analysis/subclustering/vascular/scDRS/CAD.gs \
    --zscore-file GCST90132314_gene.tsv \
    --weight zscore \
    --n-max 1000

# Using Covariate File with Study, Donor, sequencing_method
scdrs compute-score \
    --h5ad-file "/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_final.h5ad" \
    --h5ad-species human \
    --gs-file /COMMA/results/scrna/downstream_analysis/subclustering/vascular/scDRS/CAD.gs \
    --gs-species human \
    --out-folder /COMMA/results/scrna/downstream_analysis/subclustering/vascular/scDRS \
    --cov-file /COMMA/results/scrna/downstream_analysis/subclustering/vascular/scDRS/vascular.cov \
    --flag-filter-data True \
    --flag-raw-count False \
    --n-ctrl 1000 \
    --flag-return-ctrl-raw-score False \
    --flag-return-ctrl-norm-score True

# Perform Downstream Analysis
scdrs perform-downstream \
    --h5ad-file "/COMMA/results/scrna/downstream_analysis/subclustering/vascular/vascular_annotated_final.h5ad" \
    --score-file /COMMA/results/scrna/downstream_analysis/subclustering/vascular/scDRS/CAD.full_score.gz \
    --out-folder /COMMA/results/scrna/downstream_analysis/subclustering/vascular/scDRS/cell_types \
    --group-analysis vascular_cell_type \
    --gene-analysis \
    --flag-filter-data True \
    --flag-raw-count False

python /COMMA/006_GWAS_Integration/003_scDRS/Python_scripts/001_scDRS.py