#!/bin/bash

cd /COMMA/results/GWAS/s-LDSC/

wget https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90132001-GCST90133000/GCST90132314/GCST90132314_buildGRCh37.tsv
wget https://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

# Liftover hg37 to hg38 Coordinates
IN=GCST90132314_buildGRCh37.tsv
CHAIN=hg19ToHg38.over.chain.gz
module load Miniforge3
eval "$(conda shell.bash hook)"
conda activate liftover
export PATH="$PATH:conda/skylake/envs/liftover/bin"
awk 'NR==1{for(i=1;i<=NF;i++)h[$i]=i; next} {c=$h["chromosome"]; if(c!~/^chr/) c="chr"c; print c,$h["base_pair_location"]-1,$h["base_pair_location"],NR}' OFS='\t' ${IN} > input.bed
liftOver input.bed ${CHAIN} lifted.bed unmapped.bed
awk 'BEGIN{FS=OFS="\t"} FNR==NR{
  m[$4]=$1"\t"$3
  next}
FNR==1{print $0,"chromosome_hg38","base_pair_location_hg38"
  next}
(FNR in m){split(m[FNR],a,"\t")
  sub(/^chr/,"",a[1])
  print $0,a[1],a[2]}' lifted.bed ${IN} > GCST90132314_buildGRCh38.tsv

conda deactivate
conda activate ldsc
awk 'BEGIN{FS=OFS="\t"} NR==1{$10=""; sub(/\t\t/,"\t"); print $0,"SNP"; next} {$10=""; sub(/\t\t/,"\t"); print $0,$22":"$23}' GCST90132314_buildGRCh38.tsv > GCST90132314_buildGRCh38.withSNP.tsv
python conda/skylake/envs/ldsc/ldsc/munge_sumstats.py --sumstats GCST90132314_buildGRCh38.withSNP.tsv --snp SNP \
  --a1 effect_allele --a2 other_allele --p p_value --N-col n --n-min 0 \
  --info-min 0 --maf-min 0 --signed-sumstats beta,0 --out GCST90132314_buildGRCh38.munged

# Download Weights for File
cd /COMMA/reference_genome/ldsc/weights
curl -L -o https://alkesgroup.broadinstitute.org/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz/
for chr in {1..22}; do
    zcat raw/weights.hm3_noMHC.${chr}.l2.ldscore.gz |
    awk 'BEGIN{OFS="\t"}
         NR==1{print; next}
         {$2=$1":"$3; print}' |
    gzip > posid/weights.posid.${chr}.l2.ldscore.gz
done

# Convert FRQ File to per Chromosome
for chr in {1..22}; do
   awk -v c="$chr" 'BEGIN{OFS="\t"} NR==1 || $1==c' \
     /COMMA/reference_genome/ldsc/EURhg38/EUR_posid_hg38.frq \
     > /COMMA/reference_genome/ldsc/EURhg38/frq_split/EUR_posid_hg38.${chr}.frq
 done