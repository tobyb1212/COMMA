#!/bin/bash

MAKE_ANNOT="conda/skylake/envs/ldsc/ldsc/continuous_make_annot_merged.py"
LDSC="conda/skylake/envs/ldsc/ldsc/ldsc.py"
BFILE="/COMMA/reference_genome/ldsc/plink_map/merged/merged_ALL_chr"
BIM="/COMMA/reference_genome/ldsc/plink_map/merged/merged_ALL_chr.bim"
WEIGHTS_DIR="/COMMA/reference_genome/ldsc/L2/symlinks/"
InputBeds=(/COMMA/results/GWAS/s-LDSC/annotations/*.bed)
SUMSTATS="/COMMA/results/GWAS/s-LDSC/GCST90132314_buildGRCh38.munged.sumstats.gz"
OUTDIR="/COMMA/results/GWAS/s-LDSC/partitioned_heritability"
FRQ_FILE="/COMMA/reference_genome/ldsc/plink_map/merged/merged_ALL_chr"
BGBed="/COMMA/results/GWAS/s-LDSC/background/background.bed"

module load Miniforge3
eval "$(conda shell.bash hook)"
conda activate ldsc
module load BEDTools

# Preprocess BED Files
for bed in "${InputBeds[@]}"; do
  base=$(basename "${bed%.*}")
  norm="${OUTDIR}/${base}.norm.bed"
  bedtools sort -i "$bed" > "$norm"
  python "$MAKE_ANNOT" --bed-file "$norm" --bimfile "$BIM" --annot-file "${OUTDIR}/${base}.annot"
  gzip -f "${OUTDIR}/${base}.annot"
  python "$LDSC" --l2 --bfile "$BFILE" --ld-wind-kb 1000 --annot "${OUTDIR}/${base}.annot.gz" --thin-annot --out "${OUTDIR}/${base}"
done

# All unique OCRs: union of cell-type peaks and nonspecific background
{
for bed in "${InputBeds[@]}" "$BGBed"; do
    cut -f1-3 "$bed"
  done
} | bedtools sort -i - | bedtools merge -i - > "${OUTDIR}/All_OCR.norm.bed"

python "$MAKE_ANNOT" --bed-file "${OUTDIR}/All_OCR.norm.bed" --bimfile "$BIM" --annot-file "${OUTDIR}/All_OCR.annot"
gzip -f "${OUTDIR}/All_OCR.annot"
python "$LDSC" --l2 --bfile "$BFILE" --ld-wind-kb 1000 --annot "${OUTDIR}/All_OCR.annot.gz" --thin-annot --out "${OUTDIR}/All_OCR"

# Five cell types + All_OCR in one joint regression
REFS_NoBG=$(printf '%s\n' "${OUTDIR}"/*.l2.ldscore.gz |
    grep -v '/All_OCR.l2.ldscore.gz$' |
    sed 's/\.l2\.ldscore\.gz$//' | paste -sd, -)
REFS="${REFS_NoBG},${OUTDIR}/All_OCR"
ORDER_FILE="${OUTDIR}/partitioned_h2_order.txt"
echo "$REFS" | tr ',' '\n' > "$ORDER_FILE"
python "$LDSC" --h2 "$SUMSTATS" --ref-ld "$REFS" --frqfile "$FRQ_FILE" \
       --overlap-annot --w-ld-chr "$WEIGHTS_DIR" --print-coefficients \
       --print-cov --print-delete-vals \
       --out "${OUTDIR}/partitioned_h2"